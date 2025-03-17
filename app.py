import matplotlib
matplotlib.use('Agg')
from flask import Flask, request, render_template, redirect, url_for, send_from_directory
import numpy as np
from initialization import create_reference_genome, index_reference_genome
from hardy_weinberg import *
from sanger import *
from coverage import *
from visualizations import *
import os
import time

app = Flask(__name__)

# initialize reference_genome and index
reference_length = 1000
kmer_length = 3
reference_genome = create_reference_genome(reference_length)
reference_index = index_reference_genome(reference_genome, kmer_length)


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
    # Default values
    read_length = 10
    num_reads = 100
    
    if request.method == "POST":
        read_length = int(request.form.get('read_length', 5))
        num_reads = int(request.form.get('num_reads', 10))
        
        # Create Reads
        reads = create_reads(reference_genome, read_length, num_reads)
        
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

@app.route('/mutations', methods=["GET", "POST"])
def assembly():
    return render_template('mutations.html')


if __name__ == '__main__':
    if not os.path.exists('static/plots'):
        os.makedirs('static/plots')
    app.run(host="0.0.0.0", port=5000, debug=True)
