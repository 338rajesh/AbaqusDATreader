
# %%
# ==========================================================
#               IMPORTS and INPUT
# ==========================================================

import os
import re
import pickle


# %%
# ===========================================================
#                 DEFINITIONS
# ===========================================================

def get_lines(fpath: str) -> list:
    w_flag = False
    flines = []
    with open(fpath, 'r') as file:
        for line in file:
            if line.strip() == "P R O B L E M   S I Z E":
                w_flag = True
            if w_flag:
                flines.append(line)
            if line.strip() == "THE ANALYSIS HAS BEEN COMPLETED":
                w_flag = False
    return flines


def get_ele_and_ip_tags(datlines: list[str]):
    ele_numbers = []
    int_pnt_numbers = []
    for a_dat_line in datlines:
        adat_line_parts = re.split("\s+", a_dat_line)
        a_ele_num = int(adat_line_parts[0])
        a_ipnt_num = int(adat_line_parts[1])
        if a_ele_num not in ele_numbers:
            ele_numbers.append(a_ele_num)
        if a_ipnt_num not in int_pnt_numbers:
            int_pnt_numbers.append(a_ipnt_num)
    return ele_numbers, int_pnt_numbers


def get_SSA_DA_ele_op(lines: list[str]) -> dict[int, dict[int, list[float]]]:
    op_bgn_idcs = [lidx+1 for (lidx, i) in enumerate(lines)
                   if i.startswith("THE FOLLOWING TABLE IS")]
    op_end_idcs = [lidx for (lidx, i) in enumerate(
        lines) if i.startswith("MAXIMUM")]
    all_ele_vals = {}
    for (op_bgn_idx, op_end_idx) in zip(op_bgn_idcs, op_end_idcs):
        # extract element data tables
        atbl_lines = lines[op_bgn_idx:op_end_idx]
        # get the output variables
        ele_op_tab_header = atbl_lines[0]
        assert ele_op_tab_header.startswith(
            "ELEMENT"), "Element output table header mismatch"
        all_ele_vals["output_var"] = [
            i for i in ele_op_tab_header.split(" ") if len(i) != 0][3:]
        #
        data_lines = atbl_lines[2:]



        ele_tags, ip_tags = get_ele_and_ip_tags(data_lines)
        assert len(data_lines) == (len(ele_tags) * len(ip_tags)
                                   * 2), "Mismatch in the number of data lines."
        #
        counter = 0
        for ael_tag in ele_tags:
            ael_ip_vals = {}
            for (ipidx, aip_tag) in enumerate(ip_tags):
                magnitudes = re.split("\s+", data_lines[counter])
                phase_angles = re.split("\s+", data_lines[counter+1])
                assert "SSD" in phase_angles, "SSD tag missing in phase angles data row."
                magnitudes = [float(i) for i in magnitudes[2:]]
                phase_angles = [float(i) for i in phase_angles[3:]]
                counter += 2
                ael_ip_vals[aip_tag] = magnitudes + phase_angles
            all_ele_vals[ael_tag] = ael_ip_vals
    return all_ele_vals


def get_SSA_DA_nodal_op():
    return


def get_SSA_DA_a_freq_incr_data(lines: list[str]) -> tuple:
    # receives a block of frequency increment data
    full_lines = [i.strip() for i in lines if len(i.strip()) != 0]
    #
    header = [ln for ln in lines if "INCREMENT NUMBER" in ln][0]
    header_parts = [i.strip()
                    for i in header.split(" ") if len(i.strip()) != 0]
    incr_num = int(header_parts[2])
    freq = float(header_parts[-1])
    ele_op_flag = "E L E M E N T   O U T P U T" in full_lines
    nodal_op_flag = "N O D E   O U T P U T" in full_lines
    print(f"on #{incr_num} frequency incrment at frequency {freq}; elment output: {ele_op_flag}, node output: {nodal_op_flag}")
    #
    #
    ele_output_data = {}
    nodal_output_data = {}
    if ele_op_flag:
        ele_op_idx = full_lines.index("E L E M E N T   O U T P U T")
        if nodal_op_flag:
            # both nodal and element output
            node_op_idx = full_lines.index("N O D E   O U T P U T")
            if ele_op_idx < node_op_idx:
                ele_op_data_slice = full_lines[ele_op_idx:node_op_idx]
                nod_op_data_slice = full_lines[node_op_idx:]
            else:
                nod_op_data_slice = full_lines[node_op_idx:ele_op_idx]
                ele_op_data_slice = full_lines[ele_op_idx:]
            ele_output_data = get_SSA_DA_ele_op(ele_op_data_slice)
        else:
            # only element output
            ele_op_data_slice = full_lines[ele_op_idx:]
            ele_output_data = get_SSA_DA_ele_op(ele_op_data_slice)
    elif nodal_op_flag:
        # only nodal output
        node_op_idx = full_lines.index("N O D E   O U T P U T")
        nod_op_data_slice = full_lines[node_op_idx:]

    return (incr_num, freq, ele_output_data, nodal_output_data)


def get_SSA_DA_step_data(lines: list[str]) -> dict:
    freq_incr_lpos = [ln_idx for (ln_idx, ln) in enumerate(lines)
                      if ("INCREMENT NUMBER" in ln) and ("AT FREQUENCY (CYCLES/TIME)" in ln)]
    fil_pos_ends = freq_incr_lpos + [len(lines)]
    filnum_range = [(fil_pos_ends[i], fil_pos_ends[i+1])
                    for i in range(len(freq_incr_lpos))]
    step_data = {}
    for (afir_bgn, afir_end) in filnum_range:
        freq_incr_num, freq, ele_opd, node_opd = get_SSA_DA_a_freq_incr_data(
            lines=lines[afir_bgn:afir_end])
        step_data[freq_incr_num] = (freq, ele_opd, node_opd,)
    return step_data


def get_IVOL_from_SA_data(lines: list[str]) -> dict[int, list[float]]:
    full_lines = [i.strip() for i in lines if len(i.strip()) != 0]
    if not "E L E M E N T   O U T P U T" in full_lines:
        raise ValueError(
            "Element output os header is missing from Static Analysis step.")
    else:
        ele_op_header_idx = full_lines.index("E L E M E N T   O U T P U T")
        full_lines = full_lines[ele_op_header_idx:]     
        #
        op_bgn_idcs = [lidx for (lidx, i) in enumerate(full_lines)
                       if i.startswith("THE FOLLOWING TABLE IS")]
        op_end_idcs = [lidx for (lidx, i) in enumerate(full_lines)
                       if i.startswith("MAXIMUM")]
        #
        ivols: dict[int, list[float]] = {}
        for (op_bgn_idx, op_end_idx) in zip(op_bgn_idcs, op_end_idcs):
            data_lines = full_lines[(op_bgn_idx+3):op_end_idx]
            ele_tags, ip_tags = get_ele_and_ip_tags(data_lines)
            num_ip_tags = len(ip_tags)
            conter = 0
            for ael_tag in ele_tags:
                ael_ivols = []
                for ipidx in range(num_ip_tags):
                    ael_ivols.append(
                        float(re.split("\s+", data_lines[conter])[2]))
                    conter += 1
                ivols[ael_tag] = ael_ivols
        return ivols


def get_step_data(lines: list[str]):
    step_type = [get_step_type(i)
                 for i in lines if i.strip().startswith("S T E P")][0]
    print(f"STEP type: {step_type}")
    if step_type == "SteadyStateAnalysis":
        step_data = get_SSA_DA_step_data(lines)
    elif step_type == "StaticAnalysis":
        step_data = get_IVOL_from_SA_data(lines)
    return {step_type: step_data}


def get_step_type(lin: str) -> str:
    lin = lin.strip()
    if lin.endswith("S T E A D Y   S T A T E   A N A L Y S I S"):
        return "SteadyStateAnalysis"
    elif lin.endswith("S T A T I C   A N A L Y S I S"):
        return "StaticAnalysis"
    else:
        raise Warning("Unable to find the step type")


def summarize_SSD_DS_dat_file(file_path: str):
    # trims unnecessary lines at start and end
    lines: list[str] = get_lines(file_path)
    #
    step_line_pos = [ln_idx for (ln_idx, ln) in enumerate(lines)
                     if ln.strip().startswith("S T E P")]
    step_ends = step_line_pos + [len(lines)]
    step_linum_range = [(step_ends[i], step_ends[i+1])
                        for i in range(len(step_line_pos))]

    data = {}
    for (slr_bgn, slr_end) in step_linum_range:
        print("\n\n"+f"Reading step data from line {slr_bgn} to {slr_end-1}")
        data.update(get_step_data(lines=lines[slr_bgn:slr_end]))
    return data


# %%

# ==========================================================
#                       MAIN
# ==========================================================

# dat_files_dir = r"C:\Users\admin\WorkOuts\AIAA2023\PostProcess"
dat_files_dir = r"D:\Rajesh_Nakka\workdir\AIAA2023\AnalysisReadDataCycle2"
# dat_files_dir = r"C:\Users\338ra\scrap\Validation_ViscoElasticBasicCases"
results_dir = r"D:\Rajesh_Nakka\workdir\AIAA2023\AnalysisCycle2Pickles"

for a_file_name in os.listdir(dat_files_dir):
    file_extn = a_file_name.split(".")[-1]
    if file_extn.upper() == "DAT":
        file_ID = a_file_name.split(".")[0]
        DATa = summarize_SSD_DS_dat_file(
            file_path=os.path.join(dat_files_dir, a_file_name))

        with open(os.path.join(results_dir, f"SSD_DS_all_data_{file_ID}.pkl"), 'wb') as f:
            pickle.dump(DATa, f)
