using PyPlot, LsqFit

const P_VALS = 10000
const COLS_IN_FILE = 4

const L = 100:100:1000
N = L.^2
const log_L = log(L)

values = Array{Float64,2}(P_VALS, COLS_IN_FILE);
p      = Array{Float64,1}(P_VALS)
p_inf  = Array{Float64,2}(length(L),P_VALS)
s      = Array{Float64,2}(length(L),P_VALS)
s_max  = Array{Float64,1}(length(L))
p_max  = Array{Float64,1}(length(L))
R2     = Array{Float64,1}(P_VALS)

# Collecting data

for i=1:length(L)
    filename = "../../percolation/src/data/measurements_" * string(N[i]) * ".txt"
    values = readdlm(filename)

    p = values[:,1]
    p_inf[i,:] = values[:,2]
    s[i,:]   = values[:,3]
    s_max[i] = maximum(s[i,:])
    p_max[i] = p[indmax(s[i,:])]
end

# Calculating correlation

# for i=1:P_VALS
#     R2[i] = (cor(log_L, log(p_inf[:,i])))^2
# end

# Saving and plotting correlation
# writedlm("data/R2.txt",[p R2])

# plt = plot(p,R2)
# xlabel(L"$p$")
# ylabel(L"$R^2$")
# savefig("../res/square_correlation2.svg")

# Calculating best fit linear line for P_inf

const P_C_INDEX = 5001 # Index with largest R^2 value (found by inspecting datafile)
const P_C = p[P_C_INDEX] # corresponding p
const log_p_inf = log(p_inf[:,P_C_INDEX])

model(x,p) = p[1] + p[2] * x

fit = curve_fit(model, log_L,log_p_inf,[0.5,0.5])

a = fit.param[1]
b = fit.param[2]

@printf "beta/nu = %f\n" -b

# Plotting best fit line for P_inf

f(x) = a + b*x
x = linspace(2,8,1000)

plt = plot(x,f(x),"b--", label="Best fit curve")
plot(log_L, log_p_inf,"go", label="Numerical data")

xlabel(L"$\ln (\xi)$")
ylabel(L"$\ln (P_{\infty})$")
legend()

savefig("../res/square_fit_P_inf.svg")
clf()

# Calculating best fit linear line for s

const log_s_max = log(s_max)

fit = curve_fit(model, log_L,log_s_max,[0.5,0.5])

a = fit.param[1]
b = fit.param[2]

@printf "gamma/nu = %f\n" b

# Plotting best fit line for s

x = linspace(2,8,1000)

plot(x,f(x),"b--", label="Best fit curve")
plot(log_L, log_s_max,"go", label="Numerical data")

xlabel(L"$\ln (\xi)$")
ylabel(L"$\ln (\max \langle s \rangle)$")
legend()

savefig("../res/square_fit_s.svg")
clf()

# Calculating best fit linear line for s

log_diff = log(abs(p_max - P_C))

fit = curve_fit(model, log_L, log_diff,[0.5,0.5])

a = fit.param[1]
b = fit.param[2]

@printf "nu = %f\n" -1/b

# Plotting best fit line for s

x = linspace(2,8,1000)

plot(x,f(x),"b--", label="Best fit curve")
plot(log_L, log_diff,"go", label="Numerical data")

xlabel(L"$\ln (\xi)$")
ylabel(L"$\ln | p_{\rm{max}} - p_c |$")
legend()

savefig("../res/square_fit_p_max.svg")
