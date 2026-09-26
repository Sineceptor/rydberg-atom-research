"""Audit odd/even open-chain behaviour without inferring a bulk phase.

The main chain scan must run first in the same explicit output directory.
This module adds dense low-energy state diagnostics, parity-separated figures,
and an independent verification of the largest sampled odd-chain error.
"""
from pathlib import Path
import csv
import json
import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from extension_model import ground_point, chain_couplings, basis_observables
from rydberg_model import hamiltonian, occupations

ROOT = Path(__file__).resolve().parent
COLOURS = {"nearest1":"#d55e00", "nearest2":"#0072b2", "all":"#222222"}
CUTOFFS = ["nearest1", "nearest2", "all"]


def read_csv(path):
    with path.open() as file:
        rows=list(csv.DictReader(file))
    for row in rows:
        for key, value in list(row.items()):
            if key in ["cutoff", "parity"]:
                continue
            try:
                row[key]=float(value)
            except ValueError:
                pass
        row["N"]=int(row["N"])
    return rows


def write_csv(path, rows):
    with path.open("w",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def inspect_state(number, state, energy, g, delta, cutoff, level):
    """Bit-string and bond diagnostics in the full basis, not a constraint projection."""
    bits, count, old_m2, spin = basis_observables(number)
    weights=abs(state)**2
    nn11=bits[:,:-1]*bits[:,1:]
    nn00=(1-bits[:,:-1])*(1-bits[:,1:])
    n11=nn11.sum(axis=1)
    n00=nn00.sum(axis=1)
    p=weights@bits
    signs=(-1.)**np.arange(number)
    # Atom 0 is printed on the LEFT below, independently of Qiskit's labels.
    perfect_A=((np.arange(number)%2)==0).astype(int)
    perfect_B=1-perfect_A
    index_A=int(perfect_A@(2**np.arange(number)))
    index_B=int(perfect_B@(2**np.arange(number)))
    # Classical nearest-neighbour control at Omega=0, at these same g/delta.
    classical_diagonal=-delta*count+g*n11
    minimum=float(classical_diagonal.min())
    classical_ground=np.isclose(classical_diagonal,minimum,rtol=0,atol=1e-10)
    maximum_independent=(n11==0)&(count==(number+1)//2)
    product_background=4*np.sum(p*(1-p))/number**2 + (signs@(2*p-1)/number)**2
    reverse=np.array([int(format(k,f"0{number}b")[::-1],2) for k in range(2**number)])
    row={"N":number,"g":g,"delta":delta,"cutoff":cutoff,"level":level,"energy":float(energy),
         "density":float(p.mean()),"legacy_M_squared":float(weights@old_m2),
         "staggered_spin_squared":float(weights@spin**2),
         "independent_same_marginal_spin_squared":float(product_background),
         "perfect_A_probability":float(weights[index_A]),"perfect_B_probability":float(weights[index_B]),
         "perfect_pattern_total_probability":float(weights[index_A]+weights[index_B]),
         "probability_any_00_bond":float(weights@(n00>0)),"probability_any_11_bond":float(weights@(n11>0)),
         "mean_00_bond_count":float(weights@n00),"mean_11_bond_count":float(weights@n11),
         "probability_exactly_one_00_and_no_11":float(weights@((n00==1)&(n11==0))),
         "probability_both_endpoints_excited":float(weights@(bits[:,0]*bits[:,-1])),
         "classical_NN_ground_configuration_count":int(classical_ground.sum()),
         "classical_NN_ground_manifold_weight":float(weights@classical_ground),
         "maximum_independent_set_count":int(maximum_independent.sum()),
         "maximum_independent_set_weight":float(weights@maximum_independent),
         "reflection_expectation":float(np.vdot(state,state[reverse]).real)}
    profiles=[{"N":number,"g":g,"delta":delta,"cutoff":cutoff,"level":level,"atom":i,"occupation":float(p[i])} for i in range(number)]
    bonds=[{"N":number,"g":g,"delta":delta,"cutoff":cutoff,"level":level,"left_atom":i,
            "probability_00":float(weights@nn00[:,i]),"probability_11":float(weights@nn11[:,i])} for i in range(number-1)]
    top=[]
    for rank,k in enumerate(np.argsort(weights)[::-1][:10],1):
        top.append({"N":number,"g":g,"delta":delta,"cutoff":cutoff,"level":level,"rank":rank,
                    "bitstring_atom0_left":"".join(str(int(x)) for x in bits[k]),"probability":float(weights[k]),
                    "excited_atoms":int(count[k]),"count_00_bonds":int(n00[k]),"count_11_bonds":int(n11[k])})
    return row, profiles, bonds, top


def plot_parity(figdir, coarse, maxima, diagnostics, profiles, bonds):
    plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":190})
    fig,axes=plt.subplots(2,3,figsize=(13,7.2),constrained_layout=True)
    for col,g in enumerate([2.,8.,16.]):
        for parity,marker,style in [(0,"o","-"),(1,"s","--")]:
            for cutoff in ["nearest1","nearest2"]:
                select=sorted([r for r in maxima if r["g"]==g and r["N"]%2==parity and r["cutoff"]==cutoff],key=lambda r:r["N"])
                x=[r["N"] for r in select]
                label=("NN" if cutoff=="nearest1" else "NN + second")+("; even" if parity==0 else "; odd")
                axes[0,col].plot(x,[100*r["refined_sampled_max_absolute_density_error"] for r in select],marker+style,color=COLOURS[cutoff],label=label)
                axes[1,col].plot(x,[100*r["coarse_max_absolute_centre_density_error"] for r in select],marker+style,color=COLOURS[cutoff])
        axes[0,col].set(title=f"g = V₀ / Ω = {g:g}",ylabel="Mean-density error (percentage points)",xticks=range(4,13))
        axes[1,col].set(xlabel="Number of atoms N",ylabel="Interior-density error (percentage points)",xticks=range(4,13))
        for ax in axes[:,col]:ax.grid(alpha=.18)
    axes[0,0].legend(fontsize=8)
    fig.suptitle("Interaction cutoff errors depend on parity and the measured region\nUpper row: locally refined density maxima; lower row: coarse-grid interior maxima")
    fig.savefig(figdir/"05_parity_density_and_interior_errors.png",bbox_inches="tight");plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(13,4.3),constrained_layout=True)
    for ax,name,title in zip(axes,["density","staggered_spin_squared","gap"],["Mean excitation density","Centred spin diagnostic ⟨S²⟩","Lowest energy gap / Ω"]):
        for parity,marker,style in [(0,"o","-"),(1,"s","--")]:
            for cutoff in ["nearest1","all"]:
                rows=sorted([r for r in coarse if r["g"]==8 and r["delta"]==3 and r["N"]%2==parity and r["cutoff"]==cutoff],key=lambda r:r["N"])
                ax.plot([r["N"] for r in rows],[r[name] for r in rows],marker+style,color=COLOURS[cutoff],label=("NN" if cutoff=="nearest1" else "All pairs")+("; even" if parity==0 else "; odd"))
        ax.set(xlabel="Number of atoms N",ylabel=title,xticks=range(4,13));ax.grid(alpha=.18)
    axes[0].legend(fontsize=8)
    fig.suptitle("A fixed-setting parity comparison: g = 8 and Δ / Ω = 3\nThe order diagnostic is normalized equally for perfect alternating strings of either parity")
    fig.savefig(figdir/"06_fixed_setting_parity.png",bbox_inches="tight");plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(13,7),constrained_layout=True)
    for col,number in enumerate([7,8,9]):
        for level,colour in [(0,"#222222"),(1,"#a6428b")]:
            p=[r for r in profiles if r["N"]==number and r["g"]==8 and r["delta"]==3 and r["cutoff"]=="all" and r["level"]==level]
            b=[r for r in bonds if r["N"]==number and r["g"]==8 and r["delta"]==3 and r["cutoff"]=="all" and r["level"]==level]
            axes[0,col].plot([r["atom"] for r in p],[r["occupation"] for r in p],"o-",color=colour,label="Ground state" if level==0 else "First excited state")
            axes[1,col].plot([r["left_atom"] for r in b],[r["probability_00"] for r in b],"o-",color=colour)
        axes[0,col].set(title=f"N = {number}",xlabel="Atom index",ylabel="Excitation probability",ylim=(-.02,1.02))
        axes[1,col].set(xlabel="Left atom of neighbouring bond",ylabel="Probability of 00 on that bond",ylim=(-.02,.75))
        for ax in axes[:,col]:ax.grid(alpha=.18)
    axes[0,1].legend(fontsize=9,loc="upper center")
    fig.suptitle("Low-energy spatial structure at g = 8, Δ / Ω = 3, all pairs\nThe finite-chain gap separates these particular states, not necessarily two alternating patterns")
    fig.savefig(figdir/"07_low_energy_profiles_and_bonds.png",bbox_inches="tight");plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(13,4.2),constrained_layout=True)
    for ax,number in zip(axes,[7,8,9]):
        rows=[r for r in diagnostics if r["N"]==number and r["g"]==8 and r["delta"]==3 and r["cutoff"]=="all"]
        labels=["Perfect alternating\nstrings (sum)","Exactly one 00\nand no 11", "Any 11 bond"]
        for k,row in enumerate(rows):
            ax.bar(np.arange(3)+(k-.5)*.35,[row["perfect_pattern_total_probability"],row["probability_exactly_one_00_and_no_11"],row["probability_any_11_bond"]],width=.35,label="Ground" if row["level"]==0 else "First excited",color="#222222" if row["level"]==0 else "#a6428b")
        ax.set(title=f"N = {number}",xticks=range(3),xticklabels=labels,ylabel="Probability",ylim=(0,1.02));ax.tick_params(axis="x",labelsize=8)
    axes[0].legend(fontsize=9)
    fig.suptitle("Configuration weights: a diagnostic of boundary-sensitive low-energy states\ng = 8, Δ / Ω = 3, all pairs; categories are not a general phase classification")
    fig.savefig(figdir/"08_low_energy_configuration_weights.png",bbox_inches="tight");plt.close(fig)


def run(output_root=None):
    if output_root is None:
        from run_support import run_managed
        # This needs chain CSVs, so a direct invocation computes those first.
        def worker(folder):
            from run_extension import run as run_chains
            run_chains(folder)
            return run(folder)
        return run_managed(worker,"parity")
    out=Path(output_root);data=out/"results";fig=out/"figures"
    coarse=read_csv(data/"ground_states_coarse.csv")
    maxima=read_csv(data/"maximum_errors.csv")
    diagnostics,profiles,bonds,top,energy_rows,overlaps=[],[],[],[],[],[]
    # The g=2 cases are controls outside the simple 0<Delta<g blockade interval.
    settings=[(n,8.,3.) for n in range(4,13)]+[(8,2.,3.),(8,2.,4.25)]
    checks={"dense_vs_scan_energy_error":0.,"dense_vs_scan_gap_error":0.,"dense_vs_scan_density_error":0.,
            "ideal_pattern_S_squared_error":0.,"classical_maximum_independent_set_count_error":0.,
            "diagnostic_eigenpair_residual":0.,"independent_high_N_energy_error":0.,
            "independent_high_N_gap_error":0.,"independent_high_N_density_error":0.}
    for number,g,delta in settings:
        vectors_by_cutoff={}
        for cutoff in CUTOFFS:
            H=hamiltonian(chain_couplings(number,g,cutoff),delta=delta)
            if number<=9:
                energies,vectors=eigh(H.toarray(),subset_by_index=[0,7])
                method="dense"
            else:
                initial=np.random.default_rng(572119).normal(size=H.shape[0])
                energies,vectors=eigsh(H,k=8,which="SA",tol=1e-12,v0=initial,ncv=60)
                order=np.argsort(energies);energies,vectors=energies[order],vectors[:,order]
                method="independent sparse k=8"
            checks["diagnostic_eigenpair_residual"]=max(checks["diagnostic_eigenpair_residual"],max(np.linalg.norm(H@vectors[:,j]-energies[j]*vectors[:,j]) for j in range(8)))
            vectors_by_cutoff[cutoff]=vectors
            for level,energy in enumerate(energies):
                energy_rows.append({"N":number,"g":g,"delta":delta,"cutoff":cutoff,"level":level,"energy":float(energy),"energy_above_ground":float(energy-energies[0])})
            for level in [0,1]:
                row,p,b,t=inspect_state(number,vectors[:,level],energies[level],g,delta,cutoff,level)
                row["diagonalization_method"]=method
                diagnostics.append(row);profiles.extend(p);bonds.extend(b);top.extend(t)
            reference=next(r for r in coarse if r["N"]==number and r["g"]==g and r["delta"]==delta and r["cutoff"]==cutoff)
            prefix="dense_vs_scan" if number<=9 else "independent_high_N"
            checks[prefix+"_energy_error"]=max(checks[prefix+"_energy_error"],abs(energies[0]-reference["ground_energy"]))
            checks[prefix+"_gap_error"]=max(checks[prefix+"_gap_error"],abs(energies[1]-energies[0]-reference["gap"]))
            checks[prefix+"_density_error"]=max(checks[prefix+"_density_error"],abs(diagnostics[-2]["density"]-reference["density"]))
        for cutoff in CUTOFFS[:2]:
            for source_level in [0,1]:
                for all_level in [0,1]:
                    overlap=abs(np.vdot(vectors_by_cutoff[cutoff][:,source_level],vectors_by_cutoff["all"][:,all_level]))**2
                    overlaps.append({"N":number,"g":g,"delta":delta,"cutoff":cutoff,"cutoff_level":source_level,"all_pair_level":all_level,"squared_overlap":float(overlap)})
    # Enumerate both perfect strings and the classical NN ground manifold.
    pattern_rows=[]
    for number in range(4,13):
        bits,count,m2,spin=basis_observables(number)
        n11=(bits[:,:-1]*bits[:,1:]).sum(axis=1)
        maximum_independent=(n11==0)&(count==(number+1)//2)
        expected=1 if number%2 else number//2+1
        checks["classical_maximum_independent_set_count_error"]=max(checks["classical_maximum_independent_set_count_error"],abs(maximum_independent.sum()-expected))
        for label,values in [("A",(np.arange(number)%2==0).astype(int)),("B",(np.arange(number)%2==1).astype(int))]:
            index=int(values@(2**np.arange(number)))
            checks["ideal_pattern_S_squared_error"]=max(checks["ideal_pattern_S_squared_error"],abs(spin[index]**2-1))
            pattern_rows.append({"N":number,"pattern":label,"bitstring_atom0_left":"".join(map(str,values)),"legacy_M_squared":float(m2[index]),"staggered_spin_squared":float(spin[index]**2),"maximum_independent_set_count":int(maximum_independent.sum())})
    # Independently repeat the largest odd-N sampled error with another seed.
    oddmax=max([r for r in maxima if r["N"]%2==1 and r["cutoff"]=="nearest1"],key=lambda r:r["refined_sampled_max_absolute_density_error"])
    repeats=[]
    for cutoff in CUTOFFS:
        row,_,_=ground_point(oddmax["N"],oddmax["g"],oddmax["refined_max_delta"],cutoff,method="sparse",seed=991827)
        repeats.append(row)
    difference=abs(repeats[0]["density"]-repeats[2]["density"])
    checks["odd_maximum_independent_seed_density_error"]=abs(difference-oddmax["refined_sampled_max_absolute_density_error"])
    # Dense N=9 check at its own maximum, independently of the iterative solver.
    ninemax=max([r for r in maxima if r["N"]==9 and r["cutoff"]=="nearest1"],key=lambda r:r["refined_sampled_max_absolute_density_error"])
    dense_repeat=[]
    for cutoff in ["nearest1","all"]:
        row,_,_=ground_point(9,ninemax["g"],ninemax["refined_max_delta"],cutoff,method="dense")
        dense_repeat.append(row)
    dense_diff=abs(dense_repeat[0]["density"]-dense_repeat[1]["density"])
    checks["N9_maximum_dense_density_error"]=abs(dense_diff-ninemax["refined_sampled_max_absolute_density_error"])
    for name,value in checks.items():
        assert value<1e-8,(name,value)
    for name,rows in [("parity_low_states.csv",diagnostics),("parity_site_profiles.csv",profiles),("parity_bond_probabilities.csv",bonds),("parity_top_bitstrings.csv",top),("parity_low_energies.csv",energy_rows),("parity_state_overlaps.csv",overlaps),("parity_perfect_pattern_normalization.csv",pattern_rows),("parity_odd_maximum_independent_repeats.csv",repeats),("parity_N9_dense_maximum_repeats.csv",dense_repeat)]:
        write_csv(data/name,rows)
    write_csv(data/"parity_fixed_setting_scan.csv",[r for r in coarse if r["g"]==8 and r["delta"]==3])
    plot_parity(fig,coarse,maxima,diagnostics,profiles,bonds)
    summary={"scope":"Parity and boundary diagnostics in finite open chains; no bulk phase or universal mechanism claim.",
             "settings":[{"N":n,"g":g,"delta":d} for n,g,d in settings],
             "state_diagnostics":diagnostics,"largest_odd_density_error":oddmax,
             "independent_odd_repeat_density_difference":difference,"dense_N9_density_difference":dense_diff,
             "N9_maximum_setting":ninemax,"checks":checks,
             "order_definitions":{"legacy_M":"(2/N) sum_i (-1)^i n_i",
                 "centered_S":"(1/N) sum_i (-1)^i (2n_i-1)",
                 "relation":"S=M-(N mod 2)/N when index i starts at zero",
                 "independent_uniform_S_squared":"4p(1-p)/N + (2p-1)^2*(N mod 2)^2/N^2"},
             "bitstring_convention":"Atom 0 is on the left in parity_top_bitstrings.csv, unlike standard Qiskit labels.",
             "classical_interpretation":"For the NN-only model at Omega=0 and 0<delta<g, even N=2m has m+1 maximum independent sets with m excitations; odd N=2m+1 has one with m+1 excitations. Finite Omega and interaction tails change weights and splittings."}
    (data/"parity_summary.json").write_text(json.dumps(summary,indent=2))
    lines=["# Scientific audit: parity revision", "", "This audit tests the criticism with new calculations; it does not treat the criticism as numerical evidence.", "",
           "## What changes from the earlier report", "",
           "The original extension sampled only even open chains. Its size statements applied only to that subsequence and were inadequate as a general description of chain-size effects. The revised grid includes every integer N from 4 to 12 at the same g=2,8,16 and 81 detunings, plus local refinements.",
           "An even chain can form either perfect alternating pattern. It cannot have a perfect alternating pattern with both endpoints excited. A classical nearest-neighbour chain in the specific interval 0<delta<g has multiple minimum-energy bit strings for even N; this is not a universal statement about all parameters or the quantum ground state.", "",
           "## A consistent odd/even order diagnostic", "",
           "The earlier M=(2/N) sum (-1)^i n_i has an offset for odd N. Its two ideal alternating strings yield (1+1/N)^2 and (1-1/N)^2, so values above one do not mean greater-than-perfect order. That column is retained only for historical comparison.",
           "The new S=(1/N) sum (-1)^i (2n_i-1) gives S^2=1 for both perfect strings, for either parity. An independent product state still has a finite background: 4p(1-p)/N+(2p-1)^2*(N mod 2)^2/N^2 for uniform p. It is not a phase identifier or universal threshold.", "",
           "## Independent checks", "", "| Check | Observed error | Required below |", "|---|---:|---:|"]
    for name,value in checks.items():lines.append(f"| {name.replace('_',' ')} | {value:.6g} | 1e-8 |")
    lines += ["",f"Largest sampled odd-chain NN density error: {difference:.12g} ({100*difference:.8g} percentage points), N={oddmax['N']}, g={oddmax['g']:g}, delta={oddmax['refined_max_delta']:g}. A separate random starting vector reproduces it.",
              "Dense N=9 diagonalisation independently checks the maximum for that odd size. The low-energy diagnostic states for N=7,8,9 at g=8,delta=3 are also computed densely and matched to the main sparse scan.", "",
              "## Reading the state evidence", "",
              "results/parity_low_states.csv reports perfect-pattern weights, 00/11 bond events, maximum-independent-set weight and order. parity_top_bitstrings.csv identifies the dominant configurations, printed with atom 0 on the LEFT. Site/bond profiles and low-energy levels show what differs between the ground and first excited states.",
              "For even chains at the examined g=8,delta=3 setting, a small gap should not be described automatically as tunnelling between only two perfect alternating configurations. The actual weights and profiles need to be consulted; they expose the relevance of configurations with a 00 bond and excited endpoints. This supports a boundary-sensitive low-energy interpretation at that setting, not a proof that every gap minimum has the same origin.",
              "The reference background 1/N assumes independent sites at uniform occupation p=1/2. The output also gives an independent-product background with the actual site marginals; comparison with it is more specific than describing the order as merely close to 1/N. Neither comparison alone assigns a phase.",
              "The g=2 diagnostic points delta=3 and 4.25 lie outside 0<delta<g. Their 11-bond probabilities and densities test a different occupancy regime. Do not transfer the simple maximum-independent-set picture or a strong-blockade interpretation throughout that scan.",
              "Interior density excludes only the two endpoints; it does not remove all boundary influence. Parity-separated errors provide evidence of boundary sensitivity, not a complete decomposition of error into a boundary part and a bulk part.", "",
              "## Portable regression and unchanged originals", "",
              "The old N=8 CSV is copied into reference_data with SHA-256 metadata. Every main chain run checks its digest and recomputes all 62 historical points; there is no external personal path or silent skip. This guards the unchanged baseline convention. Deliberately changing that convention requires a reviewed new reference, not bypassing a failing check.",
              "Earlier even-only results are superseded for size interpretation by this corrected odd/even analysis.", ""]
    (out/"SCIENTIFIC_AUDIT.md").write_text("\n".join(lines))
    return summary


if __name__=="__main__":
    run()
