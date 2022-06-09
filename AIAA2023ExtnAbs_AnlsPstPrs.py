
# %%
from numpy import sin, cos, deg2rad, array, save
import numpy as np
import pickle
import os
import matplotlib.pyplot as plt


# %%
# ===========================================================
#                 DEFINITIONS
# ===========================================================


def get_Cmatrix_Column(
    fpath: str,
    load_dir: str,
    ):

    with open(fpath, "rb") as fp:
        data: dict = pickle.load(fp)
    #
    # Get Integration point volumes
    ivols = data["StaticAnalysis"]
    #
    # Get steady state data
    SSD_data = data["SteadyStateAnalysis"]
    real_Cmatrix_column = {}
    img_Cmatrix_column = {}
    for (afi_num, afi_data) in SSD_data.items():
        #
        a_freq: float = afi_data[0]
        ele_output_data: dict = afi_data[1]
        ele_op_var: list = ele_output_data.pop("output_var")
        num_eop_var: int = len(ele_op_var)
        #
        sig_11_idx = ele_op_var.index("PHS11")
        sig_22_idx = ele_op_var.index("PHS22")
        sig_33_idx = ele_op_var.index("PHS33")
        sig_12_idx = ele_op_var.index("PHS12")
        eps_11_idx = ele_op_var.index("PHE11")
        eps_22_idx = ele_op_var.index("PHE22")
        eps_33_idx = ele_op_var.index("PHE33")
        eps_12_idx = ele_op_var.index("PHE12")
        #
        cum_vol: float = 0.0
        cum_sig_xx_real, cum_sig_xx_img = 0.0, 0.0
        cum_sig_yy_real, cum_sig_yy_img = 0.0, 0.0
        cum_sig_zz_real, cum_sig_zz_img = 0.0, 0.0
        cum_sig_xy_real, cum_sig_xy_img = 0.0, 0.0
        cum_eps_xx_real, cum_eps_xx_img = 0.0, 0.0
        cum_eps_yy_real, cum_eps_yy_img = 0.0, 0.0
        cum_eps_zz_real, cum_eps_zz_img = 0.0, 0.0
        cum_eps_xy_real, cum_eps_xy_img = 0.0, 0.0
        #
        for (ael_tag, ael_info) in ele_output_data.items():
            ael_ivols = ivols[ael_tag]
            for (ael_aipvol, ael_ipdata) in zip(ael_ivols, ael_info.values()):
                cum_vol += ael_aipvol
                sig_xx_mag, sig_xx_pang = ael_ipdata[sig_11_idx], deg2rad(ael_ipdata[sig_11_idx+num_eop_var])
                sig_yy_mag, sig_yy_pang = ael_ipdata[sig_22_idx], deg2rad(ael_ipdata[sig_22_idx+num_eop_var])
                sig_zz_mag, sig_zz_pang = ael_ipdata[sig_33_idx], deg2rad(ael_ipdata[sig_33_idx+num_eop_var])
                sig_xy_mag, sig_xy_pang = ael_ipdata[sig_12_idx], deg2rad(ael_ipdata[sig_12_idx+num_eop_var])
                eps_xx_mag, eps_xx_pang = ael_ipdata[eps_11_idx], deg2rad(ael_ipdata[eps_11_idx+num_eop_var])
                eps_yy_mag, eps_yy_pang = ael_ipdata[eps_22_idx], deg2rad(ael_ipdata[eps_22_idx+num_eop_var])
                eps_zz_mag, eps_zz_pang = ael_ipdata[eps_33_idx], deg2rad(ael_ipdata[eps_33_idx+num_eop_var])
                eps_xy_mag, eps_xy_pang = ael_ipdata[eps_12_idx], deg2rad(ael_ipdata[eps_12_idx+num_eop_var])
                # conversion to real and imaginary
                sig_xx_real, sig_xx_img = sig_xx_mag*cos(sig_xx_pang), sig_xx_mag*sin(sig_xx_pang)
                sig_yy_real, sig_yy_img = sig_yy_mag*cos(sig_yy_pang), sig_yy_mag*sin(sig_yy_pang)
                sig_zz_real, sig_zz_img = sig_zz_mag*cos(sig_zz_pang), sig_zz_mag*sin(sig_zz_pang)
                sig_xy_real, sig_xy_img = sig_xy_mag*cos(sig_xy_pang), sig_xy_mag*sin(sig_xy_pang)
                eps_xx_real, eps_xx_img = eps_xx_mag*cos(eps_xx_pang), eps_xx_mag*sin(eps_xx_pang)
                eps_yy_real, eps_yy_img = eps_yy_mag*cos(eps_yy_pang), eps_yy_mag*sin(eps_yy_pang)
                eps_zz_real, eps_zz_img = eps_zz_mag*cos(eps_zz_pang), eps_zz_mag*sin(eps_zz_pang)
                eps_xy_real, eps_xy_img = eps_xy_mag*cos(eps_xy_pang), eps_xy_mag*sin(eps_xy_pang)
                #
                cum_sig_xx_real += (ael_aipvol*sig_xx_real)
                cum_sig_yy_real += (ael_aipvol*sig_yy_real)
                cum_sig_xy_real += (ael_aipvol*sig_xy_real)
                cum_eps_xx_real += (ael_aipvol*eps_xx_real)
                cum_eps_yy_real += (ael_aipvol*eps_yy_real)
                cum_eps_xy_real += (ael_aipvol*eps_xy_real)
                cum_sig_xx_img  += (ael_aipvol*sig_xx_img)
                cum_sig_yy_img  += (ael_aipvol*sig_yy_img)
                cum_sig_xy_img  += (ael_aipvol*sig_xy_img)
                cum_eps_xx_img  += (ael_aipvol*eps_xx_img)
                cum_eps_yy_img  += (ael_aipvol*eps_yy_img)
                cum_eps_xy_img  += (ael_aipvol*eps_xy_img)
        #
        avg_sig_xx_real, avg_sig_xx_img = cum_sig_xx_real/cum_vol, cum_sig_xx_img/cum_vol
        avg_sig_yy_real, avg_sig_yy_img = cum_sig_yy_real/cum_vol, cum_sig_yy_img/cum_vol
        avg_sig_xy_real, avg_sig_xy_img = cum_sig_xy_real/cum_vol, cum_sig_xy_img/cum_vol
        avg_eps_xx_real, avg_eps_xx_img = cum_eps_xx_real/cum_vol, cum_eps_xx_img/cum_vol
        avg_eps_yy_real, avg_eps_yy_img = cum_eps_yy_real/cum_vol, cum_eps_yy_img/cum_vol
        avg_eps_xy_real, avg_eps_xy_img = cum_eps_xy_real/cum_vol, cum_eps_xy_img/cum_vol
        #
        if load_dir == "XX":
            avg_strain_real  = avg_eps_xx_real
            # avg_strain_img  = avg_eps_xx_img            
        elif load_dir == "YY":
            avg_strain_real  = avg_eps_yy_real
            # avg_strain_img  = avg_eps_yy_img
        elif load_dir == "XY":
            avg_strain_real  = avg_eps_xy_real
            # avg_strain_img  = avg_eps_xy_img
        else:
            raise ValueError("Invalid loading direction ID")

        #
        print(f"Average strain in {load_dir} is {avg_strain_real}")
        real_Cmatrix_column[a_freq] = {
            "C1i" : avg_sig_xx_real/avg_strain_real,
            "C2i" : avg_sig_yy_real/avg_strain_real,
            "C3i" : avg_sig_xy_real/avg_strain_real,
        }
        #
        img_Cmatrix_column[a_freq] = {
            "C1i" : avg_sig_xx_img/avg_strain_real,
            "C2i" : avg_sig_yy_img/avg_strain_real,
            "C3i" : avg_sig_xy_img/avg_strain_real,
        }
    return real_Cmatrix_column, img_Cmatrix_column
    

def Cmatrices(fpaths: list[str]):
    C_real = {}
    C_img = {}
    C1i_real, C1i_img = get_Cmatrix_Column(fpaths[0], "XX")
    C2i_real, C2i_img = get_Cmatrix_Column(fpaths[1], "YY")
    C3i_real, C3i_img = get_Cmatrix_Column(fpaths[2], "XY")

    for ((a_freq1, C1i_val), (a_freq2, C2i_val), (a_freq3, C3i_val)) in zip(
        C1i_real.items(), C2i_real.items(), C3i_real.items()
    ):
        assert a_freq1 == a_freq2 == a_freq3, "Frequency mismatch"
        C_real[a_freq1] = array([
            list(C1i_val.values()),
            list(C2i_val.values()),
            list(C3i_val.values()),
        ])

    for ((a_freq1, C1i_val), (a_freq2, C2i_val), (a_freq3, C3i_val)) in zip(
        C1i_img.items(), C2i_img.items(), C3i_img.items()
    ):
        assert a_freq1 == a_freq2 == a_freq3, "Frequency mismatch"
        C_img[a_freq1] = array([
            list(C1i_val.values()),
            list(C2i_val.values()),
            list(C3i_val.values()),
        ])
    return C_real, C_img


def get_engg_properties(
    fpaths: list[str],
):
    Creal, Cimg = Cmatrices(fpaths)
    engg_prop_real = {}
    for (freq, aCRmat) in Creal.items():
        Sreal = np.linalg.inv(aCRmat)
        engg_prop_real[freq] = {
            "Exx" : 1.0/Sreal[0, 0],
            "Eyy" : 1.0/Sreal[1, 1],
            "Gxy" : 1.0/Sreal[2, 2],
            "nuxy": -1.0 * (Sreal[0, 1]/Sreal[0, 0]),
        }
    #
    engg_prop_img = {}
    for (freq, aCImat) in Cimg.items():
        Simg = np.linalg.inv(aCImat)
        engg_prop_img[freq] = {
            "Exx" : 1.0/Simg[0, 0],
            "Eyy" : 1.0/Simg[1, 1],
            "Gxy" : 1.0/Simg[2, 2],
            "nuxy": -1.0 * (Simg[0, 1]/Simg[0, 0]),
        }
    return engg_prop_real, engg_prop_img


# %%


# ==========================================================
#                       MAIN
# ==========================================================

# pkls_dir = r"C:\Users\admin\WorkOuts\AIAA2023\AnalysisPickles\C60E0"
pkls_dir = r"D:\Rajesh_Nakka\workdir\AIAA2023\PostProcessCycle2\C50E10"

fpaths = {}
print(f"Reading from directory: {pkls_dir}")
for a_fname in os.listdir(pkls_dir):
    file_path = os.path.join(pkls_dir, a_fname)
    a_fname_slices = a_fname.split("-")
    fextn = (a_fname_slices[-1]).split(".")[-1]
    real_num = int((a_fname_slices[-1]).split("_")[0][4:])
    if real_num not in fpaths.keys():
        fpaths[real_num] = ["", "", ""]    
    sim_ID = (a_fname_slices[-1]).split("_")[1]
    if sim_ID=="E22":
        fpaths[real_num][0] = file_path
    elif sim_ID=="E33":
        fpaths[real_num][1] = file_path
    elif sim_ID=="G23":
        fpaths[real_num][2] = file_path

strE22 = []
losE22 = []
strE33 = []
losE33 = []
strG23 = []
losG23 = []

for (k, v) in fpaths.items():
    print("On realization number: ", k)
    eng_prop_real, eng_prop_img = get_engg_properties(v)
    for (freq, Reng_prop) in eng_prop_real.items():
        strE22.append(Reng_prop["Exx"])
        strE33.append(Reng_prop["Eyy"])
        strG23.append(Reng_prop["Gxy"])
    for (freq, Ieng_prop) in eng_prop_img.items():
        losE22.append(Ieng_prop["Exx"])
        losE33.append(Ieng_prop["Eyy"])
        losG23.append(Ieng_prop["Gxy"])
    
strE22 = np.array(strE22)
losE22 = np.array(losE22)
strE33 = np.array(strE33)
losE33 = np.array(losE33)
strG23 = np.array(strG23)
losG23 = np.array(losG23)    


print(f"strE22, Mean: {np.mean(strE22)*1e-09} & Std: {np.std(strE22)*1e-09}")
print(f"strE33, Mean: {np.mean(strE33)*1e-09} & Std: {np.std(strE33)*1e-09}")
print(f"strG23, Mean: {np.mean(strG23)*1e-09} & Std: {np.std(strG23)*1e-09}")
print(f"loss Factor for E22, Mean: {np.mean(losE22/strE22)*1e03} & Std: {np.std(losE22/strE22)*1e03}")
print(f"loss Factor for E33, Mean: {np.mean(losE33/strE33)*1e03} & Std: {np.std(losE33/strE33)*1e03}")
print(f"loss Factor for G23, Mean: {np.mean(losG23/strG23)*1e03} & Std: {np.std(losG23/strG23)*1e03}")



# NeatMatrix_eng_prop_real, NeatMatrix_eng_prop_img = get_engg_properties([
#     r"C:\Users\admin\WorkOuts\AIAA2023\PostProcess\SSD_DS_all_data_NeatMatrixE22.pkl",
#     r"C:\Users\admin\WorkOuts\AIAA2023\PostProcess\SSD_DS_all_data_NeatMatrixE33.pkl",
#     r"C:\Users\admin\WorkOuts\AIAA2023\PostProcess\SSD_DS_all_data_NeatMatrixG23.pkl",
# ])

# RUC2D50vf_eng_prop_real, RUC2D50vf_eng_prop_img = get_engg_properties([
#     r"C:\Users\admin\WorkOuts\AIAA2023\PostProcess\SSD_DS_all_data_RUC2D_50vf_E22.pkl",
#     r"C:\Users\admin\WorkOuts\AIAA2023\PostProcess\SSD_DS_all_data_RUC2D_50vf_E33.pkl",
#     r"C:\Users\admin\WorkOuts\AIAA2023\PostProcess\SSD_DS_all_data_RUC2D_50vf_G23.pkl",
# ])


# time_const = 20.0

# # %%        Plot noramlized storage moduli E22/Em vs. (2 *pi * f * time_const)

# freqvals: list[float] = []
# normalized_stMod_E22: list[float] = []
# normalized_stMod_E33: list[float] = []
# normalized_stMod_G23: list[float] = []
# for ((f_nm, prop_nm), (f_ruc, prop_ruc)) in zip(NeatMatrix_eng_prop_real.items(), RUC2D50vf_eng_prop_real.items()):
#     assert f_nm == f_ruc, "Frequency Mismatch"
#     freqvals.append(2*np.pi*f_nm*time_const)
#     normalized_stMod_E22.append(prop_nm["Exx"]) #/prop_nm["Exx"])
#     normalized_stMod_E33.append(prop_nm["Eyy"]) #/prop_nm["Eyy"])
#     normalized_stMod_G23.append(prop_nm["Gxy"]) #/prop_nm["Gxy"])

# plt.figure(0, figsize=(10, 10))
# plt.plot(freqvals, normalized_stMod_E22 , "*" , label="E22 storage")
# plt.plot(freqvals, normalized_stMod_E33 , "*" , label="E33 storage")
# plt.xscale("log")
# plt.legend(loc="best")
# plt.savefig("ElaticStorageModuliWithFrequency")
# plt.close()
# plt.figure(0, figsize=(10, 10))
# plt.plot(freqvals, normalized_stMod_G23 , "*" , label="G23 storage")
# plt.xscale("log")
# plt.legend(loc="best")
# plt.savefig("ShearStorageModuliWithFrequency")
# plt.close()



# normalized_lsMod_E22: list[float] = []
# normalized_lsMod_E33: list[float] = []
# normalized_lsMod_G23: list[float] = []
# for ((f_nm, prop_nm), (f_ruc, prop_ruc)) in zip(NeatMatrix_eng_prop_img.items(), RUC2D50vf_eng_prop_img.items()):
#     assert f_nm == f_ruc, "Frequency Mismatch"
#     # freqvals.append(2*np.pi*f_nm*time_const)
#     normalized_lsMod_E22.append(prop_nm["Exx"]) #/prop_nm["Exx"])
#     normalized_lsMod_E33.append(prop_nm["Eyy"]) #/prop_nm["Eyy"])
#     normalized_lsMod_G23.append(prop_nm["Gxy"]) #/prop_nm["Gxy"])
    
# plt.figure(0, figsize=(10, 10))
# plt.plot(freqvals, normalized_lsMod_E22, "*" ,label="E22 loss")
# plt.plot(freqvals, normalized_lsMod_E33, "*" ,label="E33 loss")
# plt.xscale("log")
# plt.legend(loc="best")
# plt.savefig("ElaticLossModuliWithFrequency")
# plt.close()
# plt.figure(0, figsize=(10, 10))
# plt.plot(freqvals, normalized_lsMod_G23, "*" ,label="G23 loss")
# plt.xscale("log")
# plt.legend(loc="best")
# plt.savefig("ShearLossModuliWithFrequency")
# plt.close()


