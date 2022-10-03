def main():
    with open(sys.argv[1]) as fin:
        for line in fin:
            line = line.strip()
            if line.startswith("P") and '-,' in line:
                data = line.split('\t')
                path = data[2].rstrip('-').split('-,')[::-1]
                data[2] = '+,'.join(path) + "+"
                print('\t'.join(data))
            else:
                print(line)
                

if __name__ == "__main__":
    import sys
    main()