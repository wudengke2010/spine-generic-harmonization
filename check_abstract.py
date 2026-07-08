import re

with open(r'C:\Users\admin\WorkBuddy\2026-07-06-12-34-04\spine_generic_harmonization_paper.tex', 'r', encoding='utf-8') as f:
    lines = f.readlines()

# Find abstract paragraph
abstract_text = ""
in_abstract = False
for i, line in enumerate(lines):
    if '\\section*{Abstract}' in line:
        in_abstract = True
        continue
    if in_abstract:
        if '\\vspace' in line or '\\noindent' in line or '\\newpage' in line:
            break
        abstract_text += line

# Clean LaTeX for word count
clean = abstract_text.strip()
clean = re.sub(r'\\textbf\{([^}]*)\}', r'\1', clean)
clean = re.sub(r'\\[a-zA-Z]+\{([^}]*)\}', r'\1', clean)
clean = re.sub(r'\\[a-zA-Z]+', '', clean)
clean = re.sub(r'[{}$\\]', '', clean)
clean = re.sub(r'---', ' ', clean)
clean = re.sub(r'--', '-', clean)
clean = re.sub(r'\s+', ' ', clean).strip()
words = clean.split()
print(f'Abstract word count: {len(words)} words')
print(f'\nFull abstract:\n{clean}')
