storedStructure = load("largemaze_fixed.mat");
lines = storedStructure.lines;
hist = storedStructure.hist;
result = storedStructure.result;
board = storedStructure.board;
lines(2,2) = 0;

figure()
plot(lines(:,1),lines(:,2),'b', result(:,1),result(:,2),'y');

save('largemaze_fixed.mat','lines','board','hist','result')

