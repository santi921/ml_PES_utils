from copy import deepcopy

# read second line of file ./al2o3_md_select_4700.xyz"
file = open("./al2o3_md_select_4700.xyz", "r")
lines = file.readlines()
# write lines to new file 

lines_out = []
file.close()

#line_study_asap = lines[1]
print(len(lines))
for line in lines:
    if line.startswith("Lattice"):
        line_study_asap = deepcopy(line)
        ind_start = line_study_asap.find("SOAP")
        ind_end = line_study_asap.find("pbc", ind_start)
        
        line_study_asap_post = line_study_asap[:ind_start] + line_study_asap[ind_end:]
        
        # find if stress is present
        ind_stress = line_study_asap_post.find("stress")
        if ind_stress != -1:
            ind_stress_end = line_study_asap_post.find('"', ind_stress)
            ind_stress_2 = line_study_asap_post.find('"', ind_stress_end+1)
            line_study_asap_post = line_study_asap_post[:ind_stress] + line_study_asap_post[ind_stress_2:]
        
        ind_mag = line_study_asap_post.find("magmom")
        if ind_mag != -1:
            ind_mag_end = line_study_asap_post.find(' ', ind_mag)
            line_study_asap_post = line_study_asap_post[:ind_mag] + line_study_asap_post[ind_mag_end:]
        
        if not line_study_asap_post.endswith('\n'):
            line_study_asap_post += '\n'
        lines_out.append(line_study_asap_post)

        #lines_out.append(line)
    else:
        if not line.endswith('\n'):
            line += '\n'
        lines_out.append(line)
print(len(lines_out))
# write new file
file = open("./al2o3_md_select_4700_post.xyz", "w")
file.writelines(lines_out)