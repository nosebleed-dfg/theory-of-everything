import re, glob, os

path = 'C:/Users/funct/Desktop/axiom/shinies/'
cite_pat = re.compile(r'\\cite\{([^}]+)\}')
bib_pat = re.compile(r'\\bibitem\{([^}]+)\}')

for f in sorted(glob.glob(path + '*.tex')):
    with open(f, 'r', encoding='utf-8') as fh:
        content = fh.read()

    cites = cite_pat.findall(content)
    all_cited = set()
    for c in cites:
        for key in c.split(','):
            all_cited.add(key.strip())

    bibitems = set(bib_pat.findall(content))

    uncited = bibitems - all_cited
    missing = all_cited - bibitems

    basename = os.path.basename(f)
    if uncited or missing:
        print(f'{basename}:')
        if uncited:
            print(f'  Uncited bibitems: {uncited}')
        if missing:
            print(f'  Missing bibitems: {missing}')
        print()
