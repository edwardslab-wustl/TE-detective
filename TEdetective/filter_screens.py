# from TEdetective.io_functions import eprint

 #FILTER INFO
#        if total_clipped_rd >= 3 or ( (total_clipped_rd >= 1) and ( (total_clipped_rd_wpat+total_discord_rd) >= 5) ) or ( total_discord_rd >= 10 ): # L1base Sim
#        if total_clipped_rd >= 5 or ( total_clipped_rd >=3 and (total_clipped_rd_wpat + total_discord_rd >= 7)) or total_discord_rd >= 10: # and test_class_score == 4 \ CEU
#        if total_clipped_rd >= 5 or ( total_clipped_rd >=3 and (total_clipped_rd + total_discord_rd >= 10)) or total_discord_rd >= 15: # and test_class_score == 4 \ ALU
#        if total_clipped_rd >= 5 or ( total_clipped_rd >=3 and (total_clipped_rd + total_discord_rd >= 10)) or total_discord_rd >= 25: # and test_class_score == 4 \ BL6NJ
#        if total_clipped_rd >= 5 or ( total_clipped_rd >=3 and (total_clipped_rd + total_discord_rd >= 10)) or total_discord_rd >= 15: # and test_class_score == 4 \ LTR
#        if total_clipped_rd >= 3 or ( total_clipped_rd >=2 and (total_clipped_rd + total_discord_rd >= 5)) or total_discord_rd >= 10: # and test_class_score == 4 \ LTR-SUB
#        if total_clipped_rd >= 3 or ( total_clipped_rd >=2 and (total_clipped_rd + total_discord_rd >= 5)): # or total_discord_rd >= 25: # and test_class_score == 4 \ Ecat11
#            and ((total_clipped_rd + total_discord_rd)*100/cnt_rd >= rp): # and total_rd_left > 0 and total_rd_right > 0:
            
def initial_ins_filter( clipped_p, clipped_n, discord_p, discord_n, num_pat_p, num_pat_n):
    returnVal = False
    clipped_total = clipped_p + clipped_n
    discord_total = discord_p + discord_n
    num_pat_total = num_pat_n + num_pat_p
    if clipped_total >= 3:
        returnVal = True
    elif clipped_total + num_pat_total + discord_total >= 5 and clipped_total >= 1:
        returnVal = True
    elif discord_total >= 10:
        returnVal = True
    return returnVal
 
def initial_ins_filter_ceu( clipped_p, clipped_n, discord_p, discord_n, num_pat_p, num_pat_n):
    returnVal = False
    clipped_total = clipped_p + clipped_n
    discord_total = discord_p + discord_n
    num_pat_total = num_pat_n + num_pat_p
    if clipped_total >= 5:
        returnVal = True
    elif clipped_total + num_pat_total + discord_total >= 7 and clipped_total >= 3:
        returnVal = True
    elif discord_total >= 10:
        returnVal = True
    return returnVal
 
def initial_ins_filter_stringent( clipped_p, clipped_n, discord_p, discord_n, num_pat_p, num_pat_n):
    returnVal = False
    clipped_total = clipped_p + clipped_n
    discord_total = discord_p + discord_n
    num_pat_total = num_pat_n + num_pat_p
    if clipped_total >= 5:
        returnVal = True
    elif clipped_total + num_pat_total + discord_total >= 10 and clipped_total >= 3:
        returnVal = True
    elif discord_total >= 25:
        returnVal = True
    return returnVal
 
def polymorph_filter( clipped_p, clipped_n, discord_p, discord_n, num_pat_p, num_pat_n):
    returnVal = False
    clipped_total = clipped_p + clipped_n
    discord_total = discord_p + discord_n
    num_pat_total = num_pat_n + num_pat_p
    if clipped_total >= 1 and (clipped_total + discord_total) >= 3:
        returnVal = True
    elif clipped_total >= 1 and (clipped_total + num_pat_total + discord_total) >= 5:
        returnVal = True
    elif discord_total >= 5:
        returnVal = True
    return returnVal

def polymorph_filter_ceu( clipped_p, clipped_n, discord_p, discord_n, num_pat_p, num_pat_n):
    returnVal = False
    clipped_total = clipped_p + clipped_n
    discord_total = discord_p + discord_n
    num_pat_total = num_pat_n + num_pat_p
    if clipped_total >= 3:
        returnVal = True
    if clipped_total >= 1 and (clipped_total + discord_total) >= 3:
        returnVal = True
    elif clipped_total >= 1 and (clipped_total + num_pat_total + discord_total) >= 5:
        returnVal = True
    elif (clipped_total + num_pat_total + discord_total) >= 3 and discord_total >= 2:
        returnVal = True
    return returnVal
 