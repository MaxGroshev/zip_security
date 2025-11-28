import subprocess
import argparse
import os
import sys
from   subprocess import Popen, PIPE, STDOUT

# -----------------------------------------------------------------------------------------

def get_files_names(files_names, dir_path):
    # print (dir_path)
    data_files_counter = 0

    for path in os.listdir(dir_path):
        if os.path.isfile(os.path.join(dir_path, path)):
            data_files_counter += 1
            files_names.append(os.path.join(dir_path, path))

# -----------------------------------------------------------------------------------------

def run_compressor(args, compressors, data_files_names):

    for (compressor) in range(0, len(compressors)):
        extention = ".lzw"
        if (compressors[compressor] == "zstd"):
            extention = ".zstd"

        if (args.count_of_files == 1):
            source = args.source_dir + args.relative_path
            dest   = args.destination_dir + args.relative_path + extention
            pipe = Popen([compressors[compressor], "-f", "-k", source, "-o" , dest], stdout = PIPE)

            stdout_data = pipe.communicate()
            string_data = stdout_data[0].decode()
            print (string_data)

        else:
            for (i, file) in zip(range(args.count_of_files), data_files_names):

                print('————————————————————————START OF COMPRESSION———————————————————————')
                dest = file.replace(args.source_dir, args.destination_dir) + extention
                if (args.run_dict):
                    if (compressor == "zstd"):
                        pipe = Popen([compressors[compressor], "-f", "-D", args.zstd_dict ,"-k", file, "-o",  dest], stdout = PIPE)
                        stdout_data = pipe.communicate()
                    else:
                        pipe = Popen([compressors[compressor], "-f", "-D", args.my_dict ,"-k", file, "-o",  dest], stdout = PIPE)
                        stdout_data = pipe.communicate()
                else:
                    pipe = Popen([compressors[compressor], "-f", "-k", file, "-o",  dest], stdout = PIPE)
                    stdout_data = pipe.communicate()
                    string_data = stdout_data[0].decode()

                print('————————————————————————END OF COMPRESSION—————————————————————————\n\n')


# -----------------------------------------------------------------------------------------

def add_parse_arguments(parser):
    parser.add_argument("-c",  "--count_of_files", type = int, default = 19)
    parser.add_argument("-rp", "--relative_path",  type = str, default = "10.dat")
    parser.add_argument("-d",  "--run_dict",       type = bool, default = 0)
    parser.add_argument('-md',  "--my_dict", type = str, default =
                    os.path.dirname(os.path.abspath(__file__)) + '/my_dictionary/bin_dict/lzw.dict')
    parser.add_argument('-zd',  "--zstd_dict", type = str, default =
                    os.path.dirname(os.path.abspath(__file__)) + '/my_dictionary/bin_dict/zstd.dict')
    parser.add_argument('-src',"--source_dir",     type = str, default =
                    os.path.dirname(os.path.abspath(__file__)) + '/my_input_dat/bin_dat/')
    parser.add_argument('-dest',  '--destination_dir',          type = str, default =
                    os.path.dirname(os.path.abspath(__file__)) + '/my_compress_dat/')

# -----------------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    add_parse_arguments(parser)
    args = parser.parse_args()

    data_files_names = []
    compressors = ["zstd", "./lzw"]

    get_files_names(data_files_names, args.source_dir)
    run_compressor(args, compressors, data_files_names)
