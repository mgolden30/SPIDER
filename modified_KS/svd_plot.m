clf

clear

A = 2*rand(5,4)-1;

s1 = svd(A);
s2 = svd(A(1:4,:));


scatter(1:4, s1, 'filled');
hold on
%scatter(1:4, s2, 'filled');
for i = 1:numel(s2)
  yline(s2(i));
end
hold off

xlim([0.5, 4.5])
xticks([1:4]);

yticks([0:2]);
ylim([0,2]);

xlabel("n");
ylabel("\sigma_n");

saveas(gcf, 'figs/singular.png');