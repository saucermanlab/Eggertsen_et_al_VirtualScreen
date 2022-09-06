% Drug Action Figure

w = 1; n = 1.223; EC50 = 0.5; 
x = [0:0.1:1];
beta = (EC50.^n - 1)./(2*EC50.^n - 1); 
K = (beta - 1).^(1./n); 
fact= w.*(beta.*x.^n)./(K.^n + x.^n);

%% Non-competitive   

figure
subplot(2,2,1)
plot(x,fact,'k','LineWidth',2.0)
hold on
red = [1,0,0];
for i=1:5
    plot(x,fact*(1-0.2*i)+0.2*i,'Color',(2/3)^(5-i)*red+1-(2/3)^(5-i),'LineWidth',2.0)
end
hold off
xlabel('Upstream Node Activity','FontSize',10)
ylabel('Downstream Node Activity','FontSize',10)
title('Non-competitive Agonist','FontSize',12)

subplot(2,2,2)
plot(x,fact,'k','LineWidth',2.0)
hold on
blue = [0,0,1];
for i=1:5
    plot(x,fact*(1-0.2*i),'Color',(2/3)^(5-i)*blue+1-(2/3)^(5-i),'LineWidth',2.0)
end
hold off
xlabel('Upstream Node Activity','FontSize',10)
ylabel('Downstream Node Activity','FontSize',10)
title('Non-competitive Antagonist','FontSize',12)

%% Competitive 

f = zeros(length(x),5);
for d = 1:5
    f(:,d)= w.*(beta.*(x+d*0.2).^n)./(K.^n + (x+d*0.2).^n);
end
subplot(2,2,3)
plot(x,fact,'k','LineWidth',2.0)
hold on
for i=1:d
    plot(x,f(:,i),'Color',(2/3)^(5-i)*red+1-(2/3)^(5-i),'LineWidth',2.0)
end

for d = 1:5
    f(:,d)= w.*(beta.*(x-d*0.2).^n)./(K.^n + (x-d*0.2).^n);
end
ylim([0,1])
hold off
xlabel('Upstream Node Activity','FontSize',10)
ylabel('Downstream Node Activity','FontSize',10)
title('Competitive Agonist','FontSize',12)

subplot(2,2,4)
plot(x,fact,'k','LineWidth',2.0)
hold on
for i=1:d
    plot(x,f(:,i),'Color',(2/3)^(5-i)*blue+1-(2/3)^(5-i),'LineWidth',2.0)
end
ylim([0,1])
hold off
xlabel('Upstream Node Activity','FontSize',10)
ylabel('Downstream Node Activity','FontSize',10)
title('Competitive Antagonist','FontSize',12)