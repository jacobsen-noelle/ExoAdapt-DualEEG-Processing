% t_remaining() - time remaining to complete for loop
% mytic = tic; count = number of items processed so far; numFiles = t
%
function t_remaining(mytic,fin,count,total)
    % estimate time remaining
    assignin('base','count',count+1) %keep track of how many subjects we've processed
    fin(end+1) = toc(mytic);
    count = count+1;
    est = mean(fin)*(total-count)/60;
    assignin('base', 'fin',fin)
    fprintf('Estimated time remaining: %i minutes\n', round(est));
end