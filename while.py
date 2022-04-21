import sys

def write_file(file1_path, file2_path, file3_path):
    file2_dict = {}
    with open(file2_path, 'r') as file:
        for ind, line in enumerate(file):
            file2_dict[line.strip()] = ind

    output = open(file3_path, 'w')

    with open(file1_path, 'r') as file:
        for line in file:
            if line.strip().split(' ')[0] in file2_dict:
                output.write(line)

    output.close()
    return 

if "main" == __name__:
    file1_path = sys.argv[1]
    file2_path = sys.argv[2]
    file3_path = sys.argv[3]
    write_file(file1_path, file2_path, file3_path)
