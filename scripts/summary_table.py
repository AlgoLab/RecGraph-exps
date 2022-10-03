def row_with_bold(ab, rs):
    if float(ab) < float(rs):
        return f'**{ab}** | {rs}'
    elif float(rs) < float(ab):
        return f'{ab} | **{rs}**'
    else:
        return f'{ab} | {rs}'

def main(p):
    with open(p) as fin:
        lines = fin.readlines()

    genes={}
    last_gene = ''
    for line in lines:
        line = line.strip()
        if line.startswith("==>"):
            last_gene = line.split('/')[-2].rstrip('_na_aln')
            genes[last_gene] = {'abpoa':None, 'rspoa':None}
        else:
            if line.startswith('**abPOA** avg per region:'):
                d=line.split()
                genes[last_gene]["abpoa"] = d[-2]
            if line.startswith('**rsPOA** avg per region:'):
                d=line.split()
                genes[last_gene]["rspoa"] = d[-2]
    
    print('| gene | abPOA Edit (avg) | rsPOA Edit (avg) |')
    print('| --- | ---: | ---: |')
    for g in genes:
        print(f'| {g} | {row_with_bold(genes[g]["abpoa"], genes[g]["rspoa"])} |')

if __name__ == "__main__":
    import sys
    main(sys.argv[1])