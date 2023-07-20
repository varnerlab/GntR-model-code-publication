include("Include.jl")

# load the data -
data = readdlm(joinpath(_PATH_TO_SIMULATIONS, "dose_response_simulations", "results_matrix.dat"));

# first column contains the gluconate concentration
gluconate_concentration = data[:,1]
response = data[:,2:end]
μ = mean(response,dims=2)
σ = std(response,dims=2)

LB = μ .- (1.96/sqrt(1))*σ
UB = μ .+ (1.96/sqrt(1))*σ

# read dose response data
dose_response_exp = CSV.read(joinpath(_PATH_TO_DATA, "dose_response.csv"),DataFrame);
gluconate_concentration_exp = dose_response_exp[!,"Gluconate_concentration (mM)"] # mM
venus_mean_exp = dose_response_exp[!,"Mean (micromolar)"] # μM
venus_stderr_exp = dose_response_exp[!,"STD_ERR (micromolar)"] # μM




plot(log10.(gluconate_concentration),  vec(UB), color="powderblue", alpha=0.80, lw=2, label="");
plot!(log10.(gluconate_concentration),  vec(LB), color="powderblue", alpha=0.80, lw=2, label="");
plot!(log10.(gluconate_concentration),  μ, fillrange = UB, color="powderblue", alpha=0.80, label="");
plot!(log10.(gluconate_concentration),  μ, fillrange = LB, color="powderblue", alpha=0.80, label="");
plot!(log10.(gluconate_concentration),  μ, color="black", lw=2, ls=:dash, label = "mean");
scatter!(log10.(gluconate_concentration_exp), venus_mean_exp, yerr=venus_stderr_exp, label="μ ± 1.96SE", msc=:black, c=:white);
ylims!(0.2, 2.0);
xlabel!("log [Gluconate](mM)", fontsize=18)
ylabel!("[Venus] (μM)", fontsize=18)
savefig(joinpath(_PATH_TO_FIGS, "dose_response_plot.pdf"))