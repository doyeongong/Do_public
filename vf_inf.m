%%% HW3 Q3

N_a = 100;
a_up = 20.;
a_low = 0;
agrid = linspace(a_low, a_up, N_a);
beta = .96;

Pi = [.9 .1 ;
.1 .9];
Z = [.1, 1];

N_z = length(Z);
VF = zeros(N_a, N_z);
ga = zeros(N_a, N_z, 'int64');
U = zeros(N_a, N_a, N_z);
objgrid = zeros(N_a, 1);

% Fill U
for a_nx = 1:N_a
  for ap_nx = 1:N_a
    for z_nx = 1:N_z
      c = max(0, Z(z_nx) + 1.03*agrid(a_nx) - agrid(ap_nx));
      U(a_nx, ap_nx, z_nx) = log(c);
    end
  end
end

% Do VFI
err = 1.
while err > 1e-5
  VF_old = VF; % this makes a copy
  EVF = zeros(size(VF));
  for a_nx = 1:N_a
    for z_nx = 1:N_z
      for zp_nx = 1:N_z
        EVF(a_nx, z_nx) = EVF(a_nx, z_nx) + Pi(z_nx, zp_nx) * VF(a_nx, zp_nx);
      end
    end
  end
  for a_nx = 1:N_a
    for z_nx = 1:N_z
      for ap_nx = 1:N_a
        Ev = EVF(ap_nx, z_nx);
        objgrid(ap_nx) = U(a_nx, ap_nx, z_nx) + beta*Ev;
      end
      [Vs, V_nx] = max(objgrid);
      VF(a_nx, z_nx) = Vs;
      ga(a_nx, z_nx) = V_nx;
    end
  end
  err = max(max(abs(VF_old - VF)));
end

% policy function graph
plot(agrid, ga)

% Set random seed
rng(1)

% Simulate
T = 3000;
N = 1000;
Anxt = zeros(T, N, 'int64');
Znxt = zeros(T, N, 'int64');
Pithreshold = cumsum(Pi, 2);
Anxt(1, :) = 1;
Znxt(1, :) = 1;
znxtp1 = 1;
for n = 1:N
  for t = 1:(T-1)
    r = rand();
    znxt = Znxt(t, n);
    for i = 1:N_z % Quick Markov chain drawing tool
      if r < Pithreshold(znxt, i)
        znxtp1 = i;
        break;
      end
    end
    anxt = Anxt(t, n);
    anxtp1 = ga(anxt, znxt);
    Znxt(t+1, n) = znxtp1;
    Anxt(t+1, n) = anxtp1;
  end
end
Anxt = Anxt(1000:end, :);
Znxt = Znxt(1000:end, :);

%plot(Anxt(1000:1500, 1:4))
%figure()
%plot(agrid(Anxt(1000:1500, 1:4)))

A_ans=agrid(Anxt);

% Aggregate asset
E_a=mean(mean(A_ans))
