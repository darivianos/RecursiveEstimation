sumError = 0;
for ii = 1:50
    trackErrorNorm = run(1);
    sumError = sumError + trackErrorNorm;
end
sumError = sumError/50;