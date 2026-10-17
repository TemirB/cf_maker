"""Reproducible teaching examples; archival plots are extracted without redrawing."""
from pathlib import Path
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/cf-beginners-matplotlib")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pymupdf

BASE = Path(__file__).resolve().parent
FIG = BASE / "figures"
FIG.mkdir(exist_ok=True)
NAVY, TEAL, GOLD, GREY = "#142d43", "#007f86", "#ce8429", "#788897"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 18,
    "axes.labelsize": 19, "xtick.labelsize": 17, "ytick.labelsize": 17,
    "axes.formatter.use_mathtext": True,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": GREY, "axes.labelcolor": NAVY, "text.color": NAVY,
    "xtick.color": NAVY, "ytick.color": NAVY, "legend.frameon": False,
    "savefig.bbox": "tight", "pdf.fonttype": 42})

def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf", transparent=True)
    fig.savefig(FIG / f"{name}.png", dpi=180, transparent=True)
    plt.close(fig)

hc = 0.1973269804  # GeV fm
q = np.linspace(0, .22, 500)
fig, ax = plt.subplots(figsize=(7.7, 4.4))
for r, col in [(3, GOLD), (6, TEAL)]:
    ax.plot(q, 1 + .7*np.exp(-(r*q/hc)**2), color=col, lw=3,
            label=fr"$R={r}$ фм")
ax.axhline(1, color=GREY, ls="--", lw=1)
ax.set(xlabel=r"$q=|\mathbf{q}|$, ГэВ/$c$", ylabel=r"$C(q)$",
       xlim=(0,.22), ylim=(.97,1.78))
ax.legend(loc="upper right"); ax.grid(alpha=.14)
save(fig, "gaussian_radii")

fig, ax = plt.subplots(figsize=(6.7, 4.4))
p = np.linspace(-.9,.9,500)
rho = np.exp(-p*p/(2*.28**2))/(np.sqrt(2*np.pi)*.28)
ax.plot(p,rho,color=TEAL,lw=3)
ax.fill_between(p,rho,color=TEAL,alpha=.1)
for value, label, color in [(.12,r"$p_{1x}$",TEAL),(.40,r"$p_{2x}$",GOLD)]:
    ax.axvline(value,color=color,ls="--",lw=2)
    ax.text(value+.02,1.26,label,color=color)
ax.annotate("",(.12,.8),(.40,.8),arrowprops=dict(arrowstyle="<->",color=NAVY,lw=2))
ax.text(.26,.67,r"$|q_x|=0.28$",ha="center",fontsize=14)
ax.set(xlabel=r"$p_x$, ГэВ/$c$",ylabel=r"$f(p_x)$, $(\mathrm{ГэВ}/c)^{-1}$",ylim=(0,1.6))
ax.grid(alpha=.14); save(fig,"single_spectrum")

rng=np.random.default_rng(19102026)
edges=np.linspace(0,.3,31); centers=(edges[1:]+edges[:-1])/2
b_exp=centers**2*np.exp(-centers/.095)*1.7e7
c_exp=1+.65*np.exp(-(4.5*centers/hc)**2)
a=rng.poisson(b_exp*c_exp); b=rng.poisson(10*b_exp)
fig, (ax,ratio)=plt.subplots(2,1,figsize=(10,5.3),sharex=True,
    gridspec_kw={"height_ratios":[1.4,1],"hspace":.13})
ax.stairs(a,edges,color=TEAL,lw=2.4,label=r"$A$")
ax.stairs(.1*b,edges,color=GOLD,lw=2.4,label=r"$0.1 B$")
ax.set(ylabel="Пар / бин"); ax.legend(fontsize=20,loc="lower right",ncol=2)
ax.ticklabel_format(axis="y",style="sci",scilimits=(3,3));ax.grid(alpha=.13)
c=10*a/b  # N=10, since B has 10 times the pair exposure
err=c*np.sqrt(1/a+1/b)
ratio.errorbar(centers,c,yerr=err,fmt="o",ms=4,color=TEAL,capsize=2,label=r"$10 A/B$")
ratio.plot(centers,c_exp,color=NAVY,lw=1.3,alpha=.6)
ratio.axhline(1,color=GREY,ls="--");ratio.set(xlabel=r"$q=|\mathbf{q}|$, ГэВ/$c$",ylabel=r"$C(q)$",ylim=(.9,1.78))
ratio.grid(alpha=.13);save(fig,"same_mixed")
np.savetxt(FIG/"same_mixed_example.csv",np.c_[centers,a,b,c,err],delimiter=",",
    header="q_GeV_c,A,B,C_N10,poisson_error",comments="")

# Keep original bitmap labels, units, colours and legends exactly as archived.
d=pymupdf.open(BASE/"sources/balinov_2025.pdf")
for page in [8,9,13,14,15,17]:
    for j,img in enumerate(d[page-1].get_images(full=True)):
        pymupdf.Pixmap(d,img[0]).save(FIG/f"balinov_2025_p{page}_{j}.png")

# Vector crops from the STAR preprint, figure panels only, no new data points.
d=pymupdf.open(BASE/"sources/star_2005.pdf")
for name, page, rect in [("star_fig4",7,(96,115,307,301)),
                         ("star_fig14",14,(316,117,533,418))]:
    out=pymupdf.open(); clip=pymupdf.Rect(rect)
    p=out.new_page(width=clip.width,height=clip.height)
    p.show_pdf_page(p.rect,d,page,clip=clip)
    out.save(FIG/f"{name}.pdf")
print("Generated teaching plots and extracted archival figures.")
