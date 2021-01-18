cluster = gcp('nocreate');
if isempty(cluster)
    cluster = parpool([4 100]);
    cluster.IdleTimeout = 1200;
end

parfor ii = 1:48
    sprintf('Now is the %d Job', ii)
end

delete(cluster)
ii = 50;
save('TestHPCSave.mat','ii')