
#Important lines:
k=17
g=25

#Minkowski Border: z=0.0784
x3= [0,0.5]
y3= [0,-0.5]
plot(x3,y3,'k-','LineWidth', 1);

text(0.21,-0.16, 'Boundary', 'FontSize', k, 'Color', 'k','Rotation',-45);

hold on;

#Horizon: z=1,1664
x4= [0,-0.5]
y4= [0,-0.5]
plot(x4,y4,'k-','LineWidth', 1);
text(-0.35,-0.3, 'Horizon', 'FontSize', k, 'Color', 'k','Rotation',45);

hold on;

#Singularity
x = linspace(0.000784, 100, 50000);

u11 = -0.001.* ((-0.28 + sqrt(x./0.001 )).^(0.9)) .* ((0.41 + sqrt(x./0.001 )).^(1.5)) .* ((0.87 + sqrt(x./0.001 )).^(-0.44));
u22 = 0.001 .* (exp(-0.2 .* atan(3.32 .* sqrt(x./0.001) + 0.14))) .* ((-1.08 + sqrt(x./0.001)).^(0.13)) .* ((x./0.001  + 0.08 .* sqrt(x./0.001 ) + 0.09).^(0.9));
T1 = 0.5.*(tanh(u11) + tanh(u22));
R1 = 0.5.*(tanh(u11) - tanh(u22));
plot(R1, T1, '-k','LineWidth', 1);
text(-0.6,0, 'Singularity', 'FontSize', k, 'Color', 'k');

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
plot(R_values1, T_values1, 'k--','LineWidth', 1);
text(-0.13,-0.8,'Radial Mass Scale','FontSize',k,'Color','k','Rotation',90)

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

#Regions
text(-0.6,-0.2,'interior','FontSize',g,'Color','g');
hold on;
text(0.2,0.1,'Minkowski','FontSize',g,'Color','g');
hold on;
text(-0.14,-0.3,'exterior','FontSize',g,'Color','g');
hold on;

#Boundary lines
text(-0.05,0.3,'r=0','FontSize',k,'Color','k','Rotation',90);
text(0.4,-0.65,'T-R=-1','FontSize',k,'Color','k','Rotation',45);
text(-0.5,-0.55,'T+R=-1','FontSize',k,'Color','k','Rotation',-45);
xlabel('R')
ylabel('T')
title('Conformal diagram of an evaportating black hole')
grid on


hold off; % Release the current plot
