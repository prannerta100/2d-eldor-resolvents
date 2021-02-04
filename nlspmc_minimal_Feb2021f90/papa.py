import glob, os, re
files = sorted(glob.glob('*.f'))
for f in files:
    print(f)
    with open(f, 'r') as fptr:
        data = fptr.read()
    print('calls', list(set(re.findall(r'call\ (.*?)\(', data))))
    with open(f, 'r') as fptr:
        data = fptr.read()
    print('includes', list(set([re.sub(r"\'",'',x) for x in re.findall(r'include\ (.*?)\.inc', data)])))
    print('\n')

