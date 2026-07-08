import re

with open(r'C:\Users\admin\WorkBuddy\2026-07-06-12-34-04\spine_generic_harmonization_paper.tex', 'r', encoding='utf-8') as f:
    content = f.read()

# Count figure environments
fig_envs = re.findall(r'\\begin\{figure\}', content)
table_envs = re.findall(r'\\begin\{table\}', content)
print(f'Figure environments: {len(fig_envs)}')
print(f'Table environments: {len(table_envs)}')
print(f'Total: {len(fig_envs) + len(table_envs)}')

# List all figure labels
labels = re.findall(r'\\label\{(fig[^}]*)\}', content)
print(f'\nFigure labels: {labels}')

# List all table labels
tlabels = re.findall(r'\\label\{(tab[^}]*)\}', content)
print(f'Table labels: {tlabels}')

# Check for any remaining old labels
old_labels = [l for l in labels if 'study_design_a' in l or 'study_design_b' in l]
print(f'\nOld labels remaining: {old_labels}')

# Abstract word count
abs_match = re.search(r'\\section\*\{Abstract\}(.*?)(?:\\vspace)', content, re.DOTALL)
if abs_match:
    abstract = abs_match.group(1).strip()
    clean = re.sub(r'\\textbf\{([^}]*)\}', r'\1', abstract)
    clean = re.sub(r'\\[a-zA-Z]+\{([^}]*)\}', r'\1', clean)
    clean = re.sub(r'\\[a-zA-Z]+', '', clean)
    clean = re.sub(r'[{}$\\]', '', clean)
    clean = re.sub(r'---', ' ', clean)
    clean = re.sub(r'--', '-', clean)
    clean = re.sub(r'\s+', ' ', clean).strip()
    words = clean.split()
    print(f'\nAbstract word count: {len(words)} words')

# Data availability check
da_match = re.search(r'\\section\*\{Data and Code Availability\}(.*?)(?:\\section)', content, re.DOTALL)
if da_match:
    da_text = da_match.group(1)
    has_upon = 'upon publication' in da_text
    has_doi = 'DOI' in da_text or 'doi' in da_text
    has_github = 'github.com' in da_text
    print(f'\nData Availability:')
    print(f'  Contains "upon publication": {has_upon}')
    print(f'  Contains DOI: {has_doi}')
    print(f'  Contains GitHub link: {has_github}')
