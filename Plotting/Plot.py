from matplotlib import pyplot as plt
from root2mpl import *
from ROOT import *
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


file = TFile("../Merged/runs_32Ar_merged.root")
# fig, ax = plt.subplots(figsize = (9, 3))
# strip = "1"
# name = "Si5_"+strip
# name1 = "Si7_"+strip
# ax.set_yscale("log")

##########FULL
# H00 = file.Get("HStrip_"+name)
# H0 = file.Get("HStrip_"+name1)
# H0.Add(H00, 1)
# DisplayTH1D(H0, ax = ax, label="Single", color="black", rebin=4)

# ax.set_title(f"Lower Detector (strip {strip})")
# ax.set_xlim(1500, 6000)
# ax.set_ylim(10, 50000)
# ax.set_ylabel("Counts/keV")
# ax.set_xlabel("Energy [keV]")
# plt.text(3400, 2000, r"$\textbf{IAS}$", fontsize=13, color="green")
# plt.text(2350, 2500, r"$\mathbf{GT_2}$", fontsize=13, color="green")
# plt.text(2050, 2500, r"$\mathbf{GT_1}$", fontsize=13, color="green")

# plt.text(3050, 3000, r"${}^{\mathbf{148}}\textbf{Gd}$", fontsize=10, color="red")
# plt.text(5050, 2700, r"${}^{\mathbf{139}}\textbf{Pu}$", fontsize=10, color="red")
# plt.text(5400, 3000, r"${}^{\mathbf{241}}\textbf{Am}$", fontsize=10, color="red")
# plt.text(5700, 1700, r"${}^{\mathbf{244}}\textbf{Cu}$", fontsize=10, color="red")

# # Add the custom legend elements
# plt.annotate(r'$\alpha$ sources', 
#              color="red",
#              xy=(0.71, 0.83), 
#              xycoords='axes fraction', 
#              fontsize=12,
#              bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.85", lw=1))

# plt.annotate(r'${}^{32}\mathrm{Ar \; protons}$', 
#              color="green",
#              xy=(0.565, 0.83), 
#              xycoords='axes fraction', 
#              fontsize=12,
#              bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.85", lw=1))

# plt.legend(loc="upper right")
# plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.25)

# plt.savefig("Si5_1.png", dpi=300)

# plt.show()

##########ZOMMED

# fig, ax = plt.subplots(figsize = (5, 5))

# # H0001 = file.Get("HStrip_"+name2).Clone("H0001")
# H001 = file.Get("HStrip_"+name).Clone("H001")
# # H001.Add(H0001, 1)
# H01 = file.Get("HStrip_"+name1).Clone("H01")
# H01.Add(H001, 1)
# DisplayTH1D(H01, ax = ax, label="Single", color="black")

# # H0002 = file.Get("HStrip_coinc_"+name2).Clone("H0002")
# H002 = file.Get("HStrip_coinc_"+name).Clone("H002")
# # H002.Add(H0002, 1)
# H02 = file.Get("HStrip_coinc_"+name1).Clone("H02")
# H02.Add(H002, 1)
# DisplayTH1D(H02, ax = ax, label="Coincidence", color="r")

# # H0003 = file.Get("HStrip_no_coinc_"+name2).Clone("H0003")
# H003 = file.Get("HStrip_no_coinc_"+name).Clone("H003")
# # H003.Add(H0003, 1)
# H03 = file.Get("HStrip_no_coinc_"+name1).Clone("H03")
# H03.Add(H003, 1)
# DisplayTH1D(H03, ax = ax, label="No Coincidence", color="blue")
# ax.set_title(f"Lower Detector (strip {strip})")
# ax.set_xlim(3260, 3400)
# ax.set_ylim(0, 2500)
# ax.set_ylabel("Counts/keV")
# ax.set_xlabel("Energy [keV]")
# plt.legend(loc="upper left", fontsize = 12)
# plt.subplots_adjust(left=0.2, right=0.95, top=0.9, bottom=0.15)
# plt.savefig("Si5_1_z.png", dpi=300)
# plt.show()

fig, ax = plt.subplots(figsize = (6, 5))
ax.set_yscale("log")
DisplayTH1D(file.Get("HSiPMHigh"), ax=ax, label="SiPM High", color="black", title="High gain", xlabel="Energy [keV]", ylabel="Counts/4keV", rebin=4)

ax.set_xlim(0, 1500)
ax.set_ylim(10, 100000)
plt.subplots_adjust(left=0.23, right=0.95, top=0.9, bottom=0.25)
plt.savefig("SiPMhigh.png", dpi=300)

fig, ax = plt.subplots(figsize = (6, 5))
ax.set_yscale("log")
DisplayTH1D(file.Get("HSiPMLow"), ax=ax, label="SiPM Low", color="black", title="Low gain", xlabel="Energy [keV]", ylabel="Counts/20keV", rebin=20)

ax.set_xlim(0, 6500)
ax.set_ylim(10, 5000)
plt.subplots_adjust(left=0.23, right=0.95, top=0.9, bottom=0.25)
plt.savefig("SiPMlow.png", dpi=300)
plt.show()