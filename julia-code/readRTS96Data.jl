function readRTS96Data();
    # Load RTS96 data
    # uses MAT
    path = "julia-code/caseRTS96.mat"
    caseRTS96 = matread(path); # import MATLAB workspace

    # Connect relevant MATLAB variables to Julia
    bus_i = caseRTS96["bus_i"][:,1]; # vector of unique bus indices (73)
    bus = caseRTS96["bus"][:,1]; # vector of generator bus indices (99)
    Sb = caseRTS96["Sb"]; # base MVA
    bustype = caseRTS96["type"][:,1]

    Gp_long = caseRTS96["Pg"]; # conventional active power output, divide by Sb later
    Gq_long = caseRTS96["Qg"]; # conventional reactive power output, divide by Sb later

    Rp = caseRTS96["Wind"]./Sb; Rp = Rp[:,1]; # renewable active generation, has zero where no farm exists
    Rq = caseRTS96["Wind"].*(0.484/Sb); Rq = Rq[:,1];

    Dp = caseRTS96["Pd"]./Sb + Rp; # I previously subtracted wind from active demand
    Dp = Dp[:,1];
    Dq = caseRTS96["Qd"]./Sb; # reactive power demand vector
    Dq = Dq[:,1];

    Pmax = caseRTS96["Pmax"]./Sb; Pmax = Pmax[:,1]; # conventional gen. max active output
    Pmin = caseRTS96["Pmin"]./Sb; Pmin = Pmin[:,1]; # conventional gen. min active output

    Qmax = caseRTS96["Qmax"]./Sb; Qmax = Qmax[:,1]; # conventional gen. max reactive output
    Qmin = caseRTS96["Qmin"]./Sb; Qmin = Qmin[:,1]; # conventional gen. min reactive output

    Vm = caseRTS96["Vm"][:,1]; # should be vector of ones
    Va = caseRTS96["Va"][:,1]; # should be vector of zeros

    Plim = caseRTS96["rateA"]./Sb; Plim = Plim[:,1];
    Vceiling = caseRTS96["Vmax"][:,1]; # one voltage mag. ceiling constraint per node
    Vfloor = caseRTS96["Vmin"][:,1]; # one voltage mag. floor constraint per node

    Vceiling[bustype .== 2] = 1; # ceiling is 1 for voltage-controlled nodes
    Vfloor[bustype .== 2] = 1; # floor is 1 for voltage-controlled nodes

    Vg = caseRTS96["Vg"][:,1]; # voltage magnitudes from Jenny's solved ACOPF

    f = int(caseRTS96["fbus"][:,1]); # "from bus" ...
    t = int(caseRTS96["tbus"][:,1]); # ... "to bus"
    r = caseRTS96["r"][:,1]; # resistance, pu
    x = caseRTS96["x"][:,1]; # reactance, pu
    b = caseRTS96["b"][:,1];  # susceptance, pu

    # map bus numbers to 1:73
    for i = 1:length(f)
        if f[i] < 201
            f[i] -= 100
        elseif f[i] < 301
            f[i] -= 176;
        else
            f[i] -= 252;
        end
    end
    for i = 1:length(t)
        if t[i] < 201
            t[i] -= 100;
        elseif t[i] < 301
            t[i] -= 176;
        else
            t[i] -= 252;
        end
    end

    busidx = Int64[1];
    for i = 2:length(bus)
        if bus[i] < 201
            push!(busidx, bus[i]-100);
        elseif bus[i] < 301
            push!(busidx, bus[i]-176);
        else
            push!(busidx, bus[i]-252);
        end
    end

    # area 1: buses 1-24
    # area 2: buses 25-48
    # area 3: buses 49-73

    Gp = zeros(length(bus_i))
    Gq = zeros(length(bus_i))

    for i in unique(busidx)
        Gp[i] = sum(Gp_long.*[busidx .== i])/Sb;
        Gq[i] = sum(Gq_long.*[busidx .== i])/Sb;
    end
    # Now Gp and Gq reflect active and reactive generation at buses 1:73 consecutively.

    N = length(bus_i);
    Nr = length(find(Rp));
    Ng = length(find(Gp));

    Y = createY(f,t,r,x,b,true)
    # spy(Y)
    # title(L"Spy plot of $Y_{bus}$")

    return  Sb, f, t, r, x, b, Y, bustype, 
            Gp, Gq, Dp, Dq, Rp, Rq, 
            Pmax, Pmin, Qmax, Qmin, 
            Plim, Vg, Vceiling, Vfloor, 
            bus_i, N, Nr, Ng
end