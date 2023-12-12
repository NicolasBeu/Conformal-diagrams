

figure;

#lines of constant t exterior
for c_values = [0.1:0.1:10,10:1:100];
  i = 1:length(c_values)
    % Calculate the upper and lower bounds for x based on the condition
    upper_bound = c_values(i)./ 0.0784;
    lower_bound = c_values(i)./ 1.1664;

    % Generate x values within the specified range
    x = linspace(lower_bound, upper_bound, 1000);
    current_c = c_values(i);

    % Calculate the corresponding y-values for the current c
    u1 = -x .* ((-0.28 + sqrt(current_c ./ x)).^(0.9)) .* ((0.41 + sqrt(current_c ./ x)).^(1.5)) .* ((0.87 + sqrt(current_c ./ x)).^(-0.44));
    u2 = -x .* (exp(-0.2 .* atan(3.32 .* sqrt(current_c ./ x) + 0.14))) .* ((1.08 - sqrt(current_c ./ x)).^(0.13)) .* ((current_c ./ x + 0.08 .* sqrt(current_c ./ x) + 0.09).^(0.9));
    T = 0.5.*(tanh(u1) + tanh(u2));
    R = 0.5.*(tanh(u1) - tanh(u2));

    % Plot the function for the current c
    plot(R, T, '-b');

    hold on; % Keep the current plot active
end

#lines of constant t interior
for k_values = [0.1:0.1:10,10:1:100];
    j = 1:length(k_values)
    % Calculate the upper and lower bounds for x based on the condition
    upper_bound = k_values(j)./ 0.0784;
    lower_bound = k_values(j)./ 1.1664;

    % Generate x values within the specified range
    x = linspace(10e-2, lower_bound, 10000);
    current_c = k_values(i);

    % Calculate the corresponding y-values for the current c
    u1 = -x .* ((-0.28 + sqrt(current_c ./ x)).^(0.9)) .* ((0.41 + sqrt(current_c ./ x)).^(1.5)) .* ((0.87 + sqrt(current_c ./ x)).^(-0.44));
    u2 = x .* (exp(-0.2 .* atan(3.32 .* sqrt(current_c ./ x) + 0.14))) .* ((-1.08 + sqrt(current_c ./ x)).^(0.13)) .* ((current_c ./ x + 0.08 .* sqrt(current_c ./ x) + 0.09).^(0.9));
    T = 0.5.*(tanh(u1) + tanh(u2));
    R = 0.5.*(tanh(u1) - tanh(u2));

    % Plot the function for the current c
    plot(R, T, '-b');

    hold on; % Keep the current plot active
end

# lines of constant r interior region
for a_values = [0.01,0.1:0.1:2,2:1:20,20:10:200,200:100:1000];
    n = 1:length(a_values)
    % Calculate the upper and lower bounds for x based on the condition
    lower_bound = a_values(n).* 0.0784;
    upper_bound = a_values(n).* 1.1664;

    % Generate x values within the specified range
    x = linspace(upper_bound, 100, 100000);
    current_c = a_values(n);

    % Calculate the corresponding y-values for the current c
    u1 = -current_c .* ((-0.28 + sqrt(x./current_c )).^(0.9)) .* ((0.41 + sqrt(x./current_c )).^(1.5)) .* ((0.87 + sqrt(x./current_c )).^(-0.44));
    u2 = current_c .* (exp(-0.2 .* atan(3.32 .* sqrt(x./current_c) + 0.14))) .* ((-1.08 + sqrt(x./current_c )).^(0.13)) .* ((x./current_c  + 0.08 .* sqrt(x./current_c ) + 0.09).^(0.9));
    T = 0.5.*(tanh(u1) + tanh(u2));
    R = 0.5.*(tanh(u1) - tanh(u2));
    linewidth = 0.8
    % Plot the function for the current c
    plot(R, T, '-r','LineWidth', linewidth);

    hold on; % Keep the current plot active
end


# lines of constant r exterior
for b_values = [0.1:0.1:2,2:1:20,20:10:200,200:100:1000];
    m = 1:length(b_values)
    % Calculate the upper and lower bounds for x based on the condition
    lower_bound = b_values(m).* 0.0784;
    upper_bound = b_values(m).* 1.1664;

    % Generate x values within the specified range
    x = linspace(lower_bound, upper_bound, 10000);
    current_c = b_values(m);

    % Calculate the corresponding y-values for the current c
    u1 = -current_c .* ((-0.28 + sqrt(x./current_c )).^(0.9)) .* ((0.41 + sqrt(x./current_c )).^(1.5)) .* ((0.87 + sqrt(x./current_c )).^(-0.44));
    u2 = -current_c .* (exp(-0.2 .* atan(3.32 .* sqrt(x./current_c) + 0.14))) .* ((1.08 - sqrt(x./current_c )).^(0.13)) .* ((x./current_c  + 0.08 .* sqrt(x./current_c ) + 0.09).^(0.9));
    T = 0.5.*(tanh(u1) + tanh(u2));
    R = 0.5.*(tanh(u1) - tanh(u2));

    % Plot the function for the current c
    plot(R, T, '-r');

    hold on; % Keep the current plot active
end
# Minkowski lines of constant t oben
for l_values = [0.1:0.1:2];
  s = 1:length(l_values)
    x=linspace(0.01,100,10000);
    current_l = l_values(s);
     % Calculate the corresponding y-values for the current c

    u1 = x + current_l;
    u2 = -x + current_l;
    T = 0.5.*(tanh(u1)+tanh(u2));
    R = 0.5.*(tanh(u1)-tanh(u2));
    % Plot the function for the current c
    plot(R, T, '-b');

    hold on; % Keep the current plot active

end

# Minkowski lines of constant t unten
for l_values = [0.1:0.1:2];
  s = 1:length(l_values)
    upper_bound = 1.18425.*l_values(s);
    x=linspace(upper_bound,100,10000);
    current_l = -1.18425.*l_values(s);
     % Calculate the corresponding y-values for the current c

    u1 = x + current_l;
    u2 = -x + current_l;
    T = 0.5.*(tanh(u1)+tanh(u2));
    R = 0.5.*(tanh(u1)-tanh(u2));
    % Plot the function for the current c
    plot(R, T, '-b');

    hold on; % Keep the current plot active

end

# Minkowski lines of constant r
for r_values = [0.1:0.1:2,2:1:20,20:10:200];
    z = 1:length(r_values);

    lower_bound = -0.09285.*r_values(z);

    x=linspace(lower_bound,100,1000);
    current_r = 0.09285.*r_values(z);

    % Calculate the corresponding y-values for the current c

    u1 = x + current_r;
    u2 = x - current_r;
    T = 0.5.*(tanh(u1)+tanh(u2));
    R = 0.5.*(tanh(u1)-tanh(u2));
    % Plot the function for the current c
    plot(R, T, '-r');

    hold on; % Keep the current plot active
end

#Important lines:
#Horizontal blue line for t=0
x1 = [0, 1];
y1 = [0, 0];
plot(x1, y1, 'b-');
hold on;

#Minkowski Border: z=0.0784
x3= [0,0.5]
y3= [0,-0.5]
plot(x3,y3,'g-','LineWidth', 1);
hold on;

#Horizon: z=1,1664
x4= [0,-0.5]
y4= [0,-0.5]
plot(x4,y4,'g-','LineWidth', 1);
hold on;

#Singularity
x = linspace(0.000784, 100, 50000);

u11 = -0.001.* ((-0.28 + sqrt(x./0.001 )).^(0.9)) .* ((0.41 + sqrt(x./0.001 )).^(1.5)) .* ((0.87 + sqrt(x./0.001 )).^(-0.44));
u22 = 0.001 .* (exp(-0.2 .* atan(3.32 .* sqrt(x./0.001) + 0.14))) .* ((-1.08 + sqrt(x./0.001)).^(0.13)) .* ((x./0.001  + 0.08 .* sqrt(x./0.001 ) + 0.09).^(0.9));
T1 = 0.5.*(tanh(u11) + tanh(u22));
R1 = 0.5.*(tanh(u11) - tanh(u22));
#plot(R1, T1, '-k','LineWidth', 1);
hold on;

#Radial mass scale z=0
R_values1 = [];
T_values1 = [];

for current_q = [0.01:0.01:0.1, 0.1:0.1:2, 2:1:20]
    c = current_q;

    % Calculate u1 and u2 using the given formulas
    u1 = -0.9458 * c;
    u2 = -0.6409 * c;

    % Calculate R and T for the current value of c
    T = 0.5 * (tanh(u1) + tanh(u2));
    R = 0.5 * (tanh(u1) - tanh(u2));

    % Append the results to the arrays
    R_values1 = [R_values1, R];
    T_values1 = [T_values1, T];
end
darkgreen = [0, 0.4, 0];

plot(R_values1, T_values1,'color',darkgreen, '--','LineWidth', 1);
hold on;

#Border lines
plot([-1,0],[0,-1],'k-','LineWidth', 1);
hold on;
plot([0,1],[-1,0],'k-','LineWidth', 1);
hold on;
plot([1,0],[0,1],'k-','LineWidth', 1);
hold on;
plot([0,0],[1,0],'k-','LineWidth', 1);
hold on;

axis([-1, 1, -1, 1]);
xlabel('R')
ylabel('T')
grid on


hold off; % Release the current plot
