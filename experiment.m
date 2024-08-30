startT = tic;

% do dome operations
x = 124411412;
for i = 1:100000000
    x = x + 1;
end

endT = toc(startT);
disp(['Time elapsed: ' num2str(endT) ' seconds.']);

logbook = struct('gaps', [], 'vals', [], 'iteration', 12, 'status', 'optimal', 'solution', rand(10,1), 'val', 0.123, 'time', 2.2);
alg_params = struct('MaxIter', 100, 'eps', 1e-6, 'verbose', true);
problem_params = struct('actv', 0.6, 'rank', 1.1, 'ecc', 0.9, 'seed', 2, 'umin', 5, 'umax', 10, 'lmin', 2, 'lmax', 4);
log = struct('alg', 'fw_step_away', 'log', logbook, 'alg_params', alg_params, 'problem_params', problem_params);
file_name = strcat(log.alg, '_', string(log.problem_params.actv), '_', string(log.problem_params.rank), '_', string(log.problem_params.ecc), '_', string(log.problem_params.seed), '.json');

log_json = jsonencode(log, PrettyPrint=true);
fid = fopen(file_name,'w');
fprintf(fid,'%s',log_json);
fclose(fid);
