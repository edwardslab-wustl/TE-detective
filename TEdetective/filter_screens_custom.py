
def initial_ins_filter_custom( clipped_p, clipped_n, discord_p, discord_n, num_pat_p, num_pat_n):
    returnVal = False #return false = does not satisy criteria, i.e. will be filtered
    clipped_total = clipped_p + clipped_n # total clipped reads, p = plus strand, n = minus strand
    discord_total = discord_p + discord_n # total discordant reads, p = plus strand, n = minus strand
    num_pat_total = num_pat_p + num_pat_n # total polyA/T reads, p = plus strand, n = minus strand
    #change filters/logic below here as you like
    if clipped_total >= 3:
        returnVal = True
    elif clipped_total + num_pat_total + discord_total >= 5 and clipped_total >= 1:
        returnVal = True
    elif discord_total >= 10:
        returnVal = True
    return returnVal

def polymorph_filter_custom( clipped_p, clipped_n, discord_p, discord_n, num_pat_p, num_pat_n):
    returnVal = False #return false = do not filter
    clipped_total = clipped_p + clipped_n # total clipped reads, p = plus strand, n = minus strand
    discord_total = discord_p + discord_n # total discordant reads, p = plus strand, n = minus strand
    num_pat_total = num_pat_p + num_pat_n # total polyA/T reads, p = plus strand, n = minus strand
    #change filters/logic below here as you like
    if clipped_total >= 3:
        returnVal = True
    if clipped_total >= 1 and (clipped_total + discord_total) >= 3:
        returnVal = True
    elif clipped_total >= 1 and (clipped_total + num_pat_total + discord_total) >= 5:
        returnVal = True
    elif (clipped_total + num_pat_total + discord_total) >= 3 and discord_total >= 2:
        returnVal = True
    return returnVal