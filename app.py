import matplotlib
matplotlib.use('Agg')
from flask import Flask, request, render_template, redirect, url_for, send_from_directory
import numpy as np
from initialization import create_reference_genome, index_reference_genome
from hardy_weinberg import *
from sanger import *
from coverage import *
from visualizations import *
from phasing import *
from ChIP_seq import *
import os
import time

app = Flask(__name__)

@app.route('/', methods=["GET", "POST"])
def index():
    return render_template('index.html')

@app.route('/hardyweinberg', methods=["GET", "POST"])
def hardy_weinberg():
    pop_size = 1000
    
    if request.method == "POST":
        p = float(request.form.get('p_value'))
        
        theoretical = calculate_theoretical_genotypes(p)
        population_emoji, observed = calculate_observed_genotypes(p, pop_size)
        
        return render_template(
            'hardy_weinberg.html',
            p_value=p,
            population_emoji=population_emoji,
            theoretical=theoretical,
            observed=observed)
    
    else:
        return render_template('hardy_weinberg.html')
        

@app.route('/sanger', methods=["GET", "POST"])
def sanger():
    seq_len = 100
    
    if request.method == "POST":
        dd_ratio = float(request.form.get('dd_ratio', 0.1))
        num_reactions = int(request.form.get('num_reactions', 1000))
        
        # Simulate Sanger
        fragment_counts, termination_sites = simulate_sanger(seq_len, dd_ratio, num_reactions)
        
        # Statistics
        mean_length = int(np.mean(termination_sites))
        std_dev = np.std(termination_sites)
        
        # Plot
        plot_object = plot_fragment_counts(fragment_counts, dd_ratio, seq_len)
        sanger_plot_png = save_plot_to_png(plot_object, 'sanger_plot_png')
        
        # Gell
        plot_object = plot_gel_electrophoresis(fragment_counts, seq_len, dd_ratio)
        gel_png = save_plot_to_png(plot_object, 'gel_png')
        
        return render_template(
            'sanger.html', 
            dd_ratio=dd_ratio,
            num_reactions=num_reactions,
            plot_png=sanger_plot_png, 
            mean_length=mean_length, 
            std_dev=round(std_dev, 2),
            gel_png=gel_png,
            )
    else:
        return render_template(
            'sanger.html', 
            log_dd_ratio=np.log10(0.1), 
            log_num_reactions=np.log10(1000)
        )

@app.route('/coverage', methods=["GET", "POST"])
def coverage():
    # initialize reference_genome and index
    reference_length = 1000
    kmer_length = 3
    reference_genome = create_reference_genome(reference_length)
    reference_index = index_reference_genome(reference_genome, kmer_length)
    
    # Default values
    read_length = 10
    num_reads = 100
        
    if request.method == "POST":
        read_length = int(request.form.get('read_length', 5))
        num_reads = int(request.form.get('num_reads', 10))
        
        # Create Reads
        reads = create_reads(reference_genome, read_length, num_reads)
        print(reads.shape)
        
        # Align Reads
        read_starts = align_reads(reference_genome, reference_index, kmer_length, reads)
        print(min(read_starts), max(read_starts))
        # create scaffold
        scaffold = create_scaffold(reference_length, reads, read_starts)
        
        # Calculate Coverage, Unread Bases
        coverage = calculate_coverage(read_length, num_reads, reference_length)
        unread_bases = count_unread_bases(reference_genome, scaffold)
        
        # Plot
        plot_object = plot_reads(read_length=read_length, read_starts=read_starts, scaffold=scaffold, reference_length=reference_length)
        coverage_plot_png = save_plot_to_png(plot_object, 'coverage_plot_png')
        
        colored_reference_genome = color_sequence(reference_genome)
        colored_scaffold = color_sequence(scaffold)
        return render_template(
            'coverage.html', 
            read_length=read_length, 
            num_reads=num_reads, 
            plot_png=coverage_plot_png, 
            coverage=coverage, 
            unread_bases=unread_bases, 
            colored_reference_genome=colored_reference_genome, 
            colored_scaffold=colored_scaffold
        )

    else:
        return render_template(
            'coverage.html', 
            read_length=read_length, 
            num_reads=num_reads
        )

@app.route('/static/plots/<filename>')
def plot_png(filename):
    return send_from_directory('static/plots', filename)

@app.route("/phasing", methods=["GET", "POST"])
def phasing():
    sequence_len = 200
    num_reads = 100
    
    if request.method == "POST":
        error_rate = float(request.form.get("error_rate", 0.1))
        

        # Generate actual sequence
        sequence = generate_sequence(sequence_len)
        
        # Calculate reads and read_values
        reads = np.array([simulate_one_read(sequence, error_rate=error_rate) for _ in range(num_reads)])
        read_values = calculate_read_values(reads)
        
        # Generate Consensus Sequence
        consensus_sequence = generate_consensus_sequence(read_values)
        
        # Record Misreads
        accumulated_misreads, first_misread_index = record_misreads(sequence, consensus_sequence)

        # Plot and Save Illumina Read
        fig = plot_Illumina_read(sequence, read_values, consensus_sequence, accumulated_misreads, error_rate)
        plot_png = save_plot_to_png(fig, "phasing_plot")

        return render_template(
            "phasing.html",
            error_rate=error_rate,
            plot_png=plot_png,
            actual_sequence="".join(sequence),
            read_sequence="".join(consensus_sequence),
            first_misread_index=first_misread_index
        )

    return render_template("phasing.html", error_rate=0.1)


@app.route('/chipseq', methods=["GET", "POST"])
def chipseq():
    genome_length = 500
    binding_site_length = 10
    num_reads = 100
    
    if request.method == "POST":
        specificity = int(request.form.get("specificity", 4))
        ideal_locations, good_locations = generate_non_overlapping_sites(
            genome_length, binding_site_length, count=4
        )

        # Simulate genome and reads
        reference_genome, ideal_site = create_reference_genome_chip(
            genome_length=genome_length,
            binding_site_length=binding_site_length,
            ideal_site_locations=ideal_locations,
            good_site_locations=good_locations,
        )

        read_map = create_reads_chip(
            reference_genome=reference_genome,
            binding_site=ideal_site,
            binding_site_length=binding_site_length,
            antibody_specificity=specificity,
        )

        plot_object = plot_read_map(read_map)
        plot_png = save_plot_to_png(plot_object, "chipseq_plot")

        return render_template(
            "chip-seq.html",
            specificity=specificity,
            plot_png=plot_png,
            ideal_locations=ideal_locations,
            good_locations=good_locations
        )

    else:
        return render_template("chip-seq.html", specificity=2)



if __name__ == '__main__':
    if not os.path.exists('static/plots'):
        os.makedirs('static/plots')
    app.run(host="0.0.0.0", port=5000, debug=True)
