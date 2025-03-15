# Bioinformatics Concepts Website

## Overview
This project is a web-based tool designed to introduce and simulate fundamental bioinformatics concepts. Users can explore various biological and computational principles through interactive simulations. The site is built using **Flask**, **Bootstrap**, and **Matplotlib**, with dynamic visualizations to aid learning.

## Features
### Completed Simulations:
- **Sanger Sequencing**: Explore how ddNTP/dNTP ratios affect sequencing fragment distributions. Includes a histogram of fragment lengths and a simulated gel electrophoresis output.
- **Coverage Simulation**: Visualize sequencing coverage based on read length, number of reads, and genome size.

### Coming Soon:
- **Hardy-Weinberg Equilibrium**: A simulation demonstrating allele frequency stability in a non-evolving population.
- **Mutation Models**: Interactive exploration of how mutations affect alignment reads.

## Technologies Used
- **Flask**: Backend framework for handling requests and rendering pages.
- **Bootstrap**: Responsive front-end styling.
- **Matplotlib**: Generates plots for simulations.
- **Jinja2**: Templating engine for dynamic content.

## Installation
To run this project locally:
1. Clone the repository:
   ```sh
   git clone https://github.com/yourusername/bioinformatics-concepts.git
   cd bioinformatics-concepts
   ```
2. Create a virtual environment (optional but recommended):
   ```sh
   python -m venv venv
   source venv/bin/activate  # On Windows use `venv\Scripts\activate`
   ```
3. Install dependencies:
   ```sh
   pip install -r requirements.txt
   ```
4. Run the Flask application:
   ```sh
   flask run
   ```
5. Open `http://127.0.0.1:5000/` in your browser.

## Folder Structure
```
project_root/
│── static/
│   ├── styles.css
│   ├── plots/  # Stores generated images
│── templates/
│   ├── layout.html
│   ├── index.html
│   ├── sanger.html
│   ├── coverage.html
│── app.py  # Main Flask application
│── requirements.txt
│── .gitignore
│── README.md
```

## Contributing
If you’d like to contribute:
- Fork the repository
- Create a feature branch (`git checkout -b feature-name`)
- Submit a pull request

## Contact
For requests or suggestions, email [goodier.r@northeastern.edu](mailto:goodier.r@northeastern.edu).

---
🚀 **More features coming soon! Stay tuned!**

