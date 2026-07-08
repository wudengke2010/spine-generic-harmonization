"""
Sync HTML and MD files with LaTeX changes:
1. Abstract: "spanning four acquisition families" -> "across four families"
2. Abstract: "All five methods reduced" -> "All methods reduced"
3. Abstract: "and rendered multivariate" -> ", rendering multivariate"
4. Abstract: "reproducible Pareto frontier" -> "Pareto frontier"
5. Abstract: "is strongly recommended" -> "is recommended"
6. Data Availability: "will be deposited in a public repository upon publication" -> new text
7. Figure references: "Figures 1a and 1b" -> "Figure 1" (where applicable)
8. Figure captions: merge Fig1a and Fig1b captions
"""
import re

files = [
    r'C:\Users\admin\WorkBuddy\2026-07-06-12-34-04\spine_generic_harmonization_paper.html',
    r'C:\Users\admin\WorkBuddy\2026-07-06-12-34-04\spine_generic_harmonization_paper.md',
]

replacements = [
    # 1. Abstract: "spanning four acquisition families" -> "across four families"
    ('spanning four acquisition families', 'across four families'),

    # 2. Abstract: "All five methods reduced" -> "All methods reduced"
    ('All five methods reduced', 'All methods reduced'),

    # 3. Abstract: "and rendered multivariate" -> ", rendering multivariate"
    (' and rendered multivariate vendor structure non-significant',
     ', rendering multivariate vendor structure non-significant'),

    # 4. Abstract: "reproducible Pareto frontier" -> "Pareto frontier"
    ('reproducible Pareto frontier', 'Pareto frontier'),

    # 5. Abstract: "is strongly recommended" -> "is recommended"
    ('Harmonization is strongly recommended', 'Harmonization is recommended'),

    # 6. Data Availability
    ('Full analysis code, raw and harmonized biomarker tables, and figure-generation scripts will be deposited in a public repository upon publication.',
     'All processing scripts, harmonization implementations, evaluation code, and figure-generation pipelines are archived at https://github.com/yilinzhu/spine-generic-harmonization (DOI: 10.5281/zenodo.XXXXXXX). Raw and harmonized biomarker tables are provided as Supplementary Table S1 and in the same repository.'),

    # 7. Figure references - update cross-references
    ('Figures 1a and 1b', 'Figure 1'),
    ('Figures~1a and 1b', 'Figure 1'),
    ('Figure 1a and 1b', 'Figure 1'),
    ('Fig. 1a and 1b', 'Fig. 1'),

    # 8. Merge figure captions - update Fig1a caption to include Fig1b content
    # This is handled differently in HTML vs MD, so we'll handle the caption merge manually
]

for fpath in files:
    print(f'\nProcessing: {fpath}')
    with open(fpath, 'r', encoding='utf-8') as f:
        content = f.read()

    original = content
    for old, new in replacements:
        count = content.count(old)
        if count > 0:
            content = content.replace(old, new)
            print(f'  Replaced ({count}x): "{old[:60]}..." -> "{new[:60]}..."')
        else:
            print(f'  Not found: "{old[:60]}..."')

    if content != original:
        with open(fpath, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f'  -> File updated')
    else:
        print(f'  -> No changes needed')

print('\nDone!')
