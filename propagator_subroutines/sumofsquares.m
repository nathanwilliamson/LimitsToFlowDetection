function val = sumofsquares(logDelta,logMSD,param)

logMSDmodel=param(1)+param(2)*logDelta;

val=sum( (logMSD-logMSDmodel).^2);

end