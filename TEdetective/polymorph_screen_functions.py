from collections import defaultdict
from TEdetective.io_functions import eprint
import TEdetective.polymorph_screen_filters as pm_filter

def write_results_mask(results, filter_cnt, filter_file_names, out_file):
    header = "\t".join(['Type','Chr','initial_guess', 'pass_initial_ins_filter'])
    for i in range(1,filter_cnt+1):
        name = str(i)
        if i in filter_file_names:
            name = filter_file_names[i]
        header += "\t" + name + '_filter'
    with open(out_file, 'w') as fh:
        fh.write(header + "\n")
        for key in results:
            fh.write(results[key] + "\n")
    return

def write_results(results, filter_cnt, file_name, out_file):
    total_pass = 0
    with open(out_file, 'w') as out_fh:
        with open(file_name, 'r') as in_fh:
            for line in in_fh:
                if line.startswith("#") or line.startswith("Type"):
                    out_fh.write(line)
                else:
                    line = line.strip()
                    line_data = line.split()
                    chrom = line_data[1]
                    ini_pos = line_data[2]
                    key = chrom + '_' + ini_pos
                    filterVal = check_filters( results[key], filter_cnt)
                    #eprint(key, filterVal)
                    if filterVal == False:
                        out_fh.write(line + "\n")
                        total_pass += 1
    return total_pass

def check_filters (filter_result, filter_cnt):
    filterVal = False
    if filter_result[0] == False:
        filterVal = True
    else:
        for i in range(1,filter_cnt+1):
            if filter_result[i] == True:
                filterVal = True
                continue
    return filterVal

def write_stats(total_initial_pass,total_pass,total_not_found,total_initial_predictions,total_filtered,filter_name_dict,file_name):
    with open(file_name, "w") as fh:
        fh.write('total_initial_predictions: ' + str(total_initial_predictions) + "\n")
        fh.write('total_not_found: ' + str(total_not_found) + "\n")
        for key in total_filtered:
            if key == 0:
                #fh.write('initial predictions passing base filter: ' + str(total_filtered[key]) + "\n")
                fh.write('initial predictions passing base filter: ' + str(total_initial_pass)+ "\n")
            else:
                name = str(key)
                if key in filter_name_dict:
                    name = filter_name_dict[key]
                fh.write('initial predictions found in ' + name + ': ' + str(total_filtered[key]) + "\n")
        fh.write('Total passing all filters: ' + str(total_pass) + "\n")
    return
            
                
def calc_filter_results(file_name, filter_cnt, filter_results):
    not_found_val = 'NA'
    total_initial_predictions = 0
    total_not_found = 0
    total_initial_pass = 0
    total_filtered = dict()
    total_pass = 0
    results=dict()
    for i in range(0,filter_cnt+1):
        total_filtered[i] = 0
        not_found_val += '\tNA'
    with open(file_name, 'r') as FH:
        line_count = 0
        for line in FH:
            line_count += 1
            if line_count > 1:
                line = line.strip()
                line_data = line.split()
                chrom = line_data[1]
                ini_pos = line_data[2]
                key = chrom + '_' + ini_pos
                total_initial_predictions += 1
                ini_filter_pass = False
                polymorph_filter = False
                out_line = "\t".join([line_data[0],chrom,ini_pos]) 
                if key in filter_results:
                    out_line += "\t" + "\t".join([str(x) for x in filter_results[key]])
                    for i,val in enumerate(filter_results[key]):
                        #eprint(key,i,val)
                        if i == 0 and val == True:
                            ini_filter_pass = True
                            total_initial_pass += 1
                        elif val == True:
                            polymorph_filter = True
                            total_filtered[i] += 1
                    if ini_filter_pass == True and polymorph_filter == False:
                        total_pass += 1
                else:
                    out_line += "\t" + not_found_val
                    total_not_found += 1
                results[key] = out_line
    return total_initial_pass,total_pass,total_not_found,total_initial_predictions,total_filtered, results
    

def add_filter_data (filter_input, file_name, file_num, qual_threshold, filter):
    filter_header,filter_clipped_n,filter_clipped_p,filter_discord_p,filter_discord_n,filter_num_pat_p,filter_num_pat_n = read_results_file(file_name, qual_threshold)
    for key in filter_input.keys():
        filterVal = 'NA'
        if key in filter_clipped_n:
            if filter == 'ceu':
                filterVal = pm_filter.polymorph_filter_ceu( filter_clipped_p[key], filter_clipped_n[key],
                                                  filter_discord_p[key], filter_discord_n[key],
                                                  filter_num_pat_p[key], filter_num_pat_n[key])
            else:
                filterVal = pm_filter.polymorph_filter( filter_clipped_p[key], filter_clipped_n[key],
                                              filter_discord_p[key], filter_discord_n[key],
                                              filter_num_pat_p[key], filter_num_pat_n[key])
        filter_input[key].append(filterVal)
        #filter_input[key].insert(file_num, filterVal)
    return filter_input
 
    
def filter_input_file (fileName, filter, qual_threshold):
    filter_input = defaultdict(list)
    input_header, input_clipped_n, input_clipped_p, input_discord_p, input_discord_n, input_num_pat_p, input_num_pat_n = read_results_file(fileName, qual_threshold)
    for key in input_clipped_n.keys():
        if filter == 'ceu':
            #eprint(fileName, ",", key)
            filterVal = pm_filter.initial_ins_filter_ceu(input_clipped_p[key], input_clipped_n[key], input_discord_p[key], input_discord_n[key], input_num_pat_p[key], input_num_pat_n[key])
        elif filter == 'stringent':
            filterVal = pm_filter.initial_ins_filter_stringent(input_clipped_p[key], input_clipped_n[key], input_discord_p[key], input_discord_n[key], input_num_pat_p[key], input_num_pat_n[key])
        else: 
            filterVal = pm_filter.initial_ins_filter(input_clipped_p[key], input_clipped_n[key], input_discord_p[key], input_discord_n[key])
        filter_input[key] = [filterVal]
    return filter_input


def read_results_file (fileName, quality_threshold):
    results_clipped_p = dict()
    results_clipped_n = dict()
    results_discord_p = dict()
    results_discord_n = dict()
    results_num_pat_n = dict()
    results_num_pat_p = dict()
    header = ''
    line_count = 0
    with open(fileName, 'r') as FH:
        for line in FH:
            line_count += 1
            if line_count == 1:
                header = line
            else:
                line_data = line.strip().split()
                chrom = line_data[1]
                ini_pos = line_data[2]
                key = chrom + '_' + ini_pos
                if (line_data[8] != 'NA' and float(line_data[8]) > quality_threshold):
                    results_clipped_p[key] = int(line_data[9])
                else:
                    results_clipped_p[key] = 0
                if (line_data[12] != 'NA' and float(line_data[12]) > quality_threshold):
                    results_clipped_n[key] = int(line_data[13])
                else:
                    results_clipped_n[key] = 0
                if (line_data[22] != 'NA' and float(line_data[22]) > quality_threshold):
                    results_discord_p[key] = int(line_data[23])
                else:
                    results_discord_p[key] = 0
                if (line_data[26] != 'NA' and float(line_data[26]) > quality_threshold):
                    results_discord_n[key] = int(line_data[27])
                else:
                    results_discord_n[key] = 0
                results_num_pat_n[key] = int(line_data[20])
                results_num_pat_p[key] = int(line_data[19])
    return header,results_clipped_p,results_clipped_n,results_discord_p,results_discord_n,results_num_pat_p,results_num_pat_n
 
def add_filter_alt_chrom (filter_input):
    count = 0
    for key in filter_input.keys():
        filterVal = False
        if len(key.split('_')) > 2:
            filterVal = True
            count += 1
        filter_input[key].append(filterVal)
    return filter_input, count

def add_filter_existing_data (filter_input, rmsk_file, file_name, te_type):
    index_size = 50000
    filter_dict,te_info = read_rmsk_file(rmsk_file, te_type,index_size)
    with open(file_name, 'r') as FH:
        line_count = 0
        for line in FH:
            filterVal = False
            line_count += 1
            if line_count > 1:
                line = line.strip()
                line_data = line.split()
                chrom = line_data[1]
                ini_pos = int(line_data[2])
                guess_pos = int(line_data[3])
                key = chrom + '_' + str(ini_pos)
                idx = int(float(ini_pos)/float(index_size))
                for index in [idx - 1, idx, idx + 1]:
                    if (chrom,index) in filter_dict:
                        for te_id in filter_dict[(chrom,index)]:
                            if filterVal == False:
                                te_start = te_info[te_id][0]
                                te_stop = te_info[te_id][1]
                                if te_start > te_stop:
                                    te_start = te_info[te_id][1]
                                    te_stop = te_info[te_id][0]
                                if (te_start <= ini_pos and ini_pos <= te_stop):
                                    filterVal = True
                            elif filterVal == True:
                                break
                    if filterVal == True:
                        break
                filter_input[key].append(filterVal)
    return filter_input

def read_rmsk_file (fileName, te_select, index_size):
    filter_dict = defaultdict(list)
    te_info = defaultdict(list)
    te_id = 0
    with open(fileName, 'r') as FH:
        for line in FH:
            if not line.startswith("#"):
                items = line.strip().split()
                te_type = items[3]
                if te_type == te_select:
                    te_id += 1
                    chrom = items[0]
                    te_start = int(items[1])
                    te_end = int(items[2])
                    index = int(float(te_start) / float(index_size))
                    te_info[te_id] = [te_start,te_end]
                    if (chrom,index) in filter_dict:
                        filter_dict[(chrom,index)].append(te_id)
                    else:
                        filter_dict[(chrom,index)] = [te_id]
    return filter_dict, te_info
