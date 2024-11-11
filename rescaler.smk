import pandas as pd
import subprocess

rule lines_to_pull_and_file_cleaning:
    input:
        file1="tmp/shuf.a.bed",
        file2="tmp/shuf.b.bed",
        file3="metadata/reference.hist"
    output:
        filel="lines_to_pull.tsv",
        file2="shuf.merged.tsv"
    run:
        temp2 = "normalized_total.tsv"

        with open(input.file1, 'r') as input_file, open(output.file2, "w") as output_file:
            for line in input_file:
                line = line.strip()
                items = line.split("\t")
                output_file.write(f"{items[0]}\t{items[1]}\t{items[2]}\t{items[3]}\t{items[5]}\n")

        with open(input.file2, 'r') as input_file, open(output.file2, "a") as output_file:
            for line in input_file:
                line = line.strip()
                items = line.split("\t")
                output_file.write(f"{items[0]}\t{items[1]}\t{items[2]}\t{items[4]}\t{items[5]}\n")

        frag_dict = {}
        with open(output.file2) as file:
            for line in file:
                line = line.strip()
                new_list = line.split("\t")
                bp_len = new_list[3]
                if bp_len in frag_dict:
                    frag_dict[bp_len] += 1
                else:
                    frag_dict[bp_len] = 1

        total_sum = sum(frag_dict.values())
        for key in frag_dict:
            normalized = frag_dict[key] / total_sum
            frag_dict[key] = [frag_dict[key], normalized]

        with open(input.file3) as file:
            for line in file:
                line = line.strip()
                contents = line.split("\t")
                key = contents[0]
                value = float(contents[1])
                if key in frag_dict:
                    count, normalized = frag_dict[key]
                    frag_dict[key] = [count, normalized, value]

        with open(temp2, 'w') as file:
            file.write("bp\tcount\tnormalized\treference\n")
            for key, (value1, value2, value3) in frag_dict.items():
                file.write(f"{key}\t{value1}\t{value2}\t{value3}\n")

        max_diff = 0
        max_query_freq = 0
        max_ref_freq = 0
        with open(temp2, "r") as file:
            next(file)  
            for line in file:
                line = line.strip()
                new_list = line.split("\t")
                query_freq = float(new_list[2])
                reference_freq = float(new_list[3])
                diff = reference_freq - query_freq
                if diff > max_diff:
                    max_diff = diff
                    max_query_freq = query_freq
                    max_ref_freq = reference_freq

        relative_freq = max_ref_freq
        df = pd.read_csv(temp2, sep="\t", header=0)
        df["lines_to_pull"] = (df["count"] * (df["reference"] / relative_freq)).astype(int)
        df[["bp", "lines_to_pull"]].to_csv(output[0], sep="\t", index=False, header=None)


rule sub_sampler:
    input:
        file1="shuf.merged.tsv",
        file2="lines_to_pull.tsv"
    output:
        file1=expand("sub_sampled_{nums}.tsv", nums=[i for i in range(1, 11)])
    run:
        count_dict = dict()
        with open(input.file2, 'r') as file:
            for line in file:
                line = line.strip()
                items = line.split("\t")
                bp = str(items[0])
                count = int(items[1])
                count_dict[bp] = count

        for file in output.file1:
            temp_file = "shuffled_shuf.merged.bed"
            subprocess.run(f"shuf {input.file1} > {temp_file}", shell=True)
            with open(temp_file, 'r') as input_file, open(file, 'w') as output_file:
                for line in input_file:
                    line = line.strip()
                    items = line.split("\t")
                    bp = str(items[3])
                    if (bp in count_dict.keys()) and (count_dict[bp] > 0):
                        output_file.write(line)
                        output_file.write("\n")
                        count_dict[bp] -= 1
            subprocess.run(f"rm {temp_file}", shell=True)

rule graph_data_generator:
    input: "sub_sampled_1.tsv"
    output: "sub_sampled_1_plot.tsv"
    run:
        frag_dict = dict()
        with open(input[0]) as file:
            for line in file:
                line = line.strip()
                new_list = line.split("\t")
                bp_len = new_list[3]
                if bp_len in frag_dict:
                    frag_dict[bp_len] += 1
                else:
                    frag_dict[bp_len] = 1

        total_sum = sum(frag_dict.values())
        for key in frag_dict:
            normalized = frag_dict[key] / total_sum
            frag_dict[key] = float(normalized)

        with open(output[0], 'w') as file:
            file.write("bp\tnormalized\n")
            for key, value in frag_dict.items():
                file.write(f"{key}\t{value}\n")
