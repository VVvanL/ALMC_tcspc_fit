function updateProgress(total, start_tic)
    persistent count prev_percent last_msg_len
    if isempty(count)
        count = 0; 
        prev_percent = -1;
        last_msg_len = 0;
    end
    
    count = count + 1;
    current_percent = floor((count/total) * 100);
    
    % Update every 1% or on the very last item
    if current_percent > prev_percent || count == total
        prev_percent = current_percent;
        
        % Calculate Elapsed and Remaining
        elapsed = toc(start_tic);
        remaining = (elapsed / (count/total)) - elapsed;

        % Build Progress Bar
        bar_width = 25;
        filled = floor((count/total) * bar_width);
        bar_str = [repmat('█', 1, filled), repmat('_', 1, bar_width - filled)];
        
        % Construct Message
        % Format: [███___] | 45% | Elapsed: 01:20 | Est. Left: 00:45
        msg = sprintf('[%s] | %d%% | Elapsed: %s | Est. Left: %s  ', ...
            bar_str, current_percent, formatTime(elapsed), formatTime(remaining));
        
        % Overwrite previous line
        fprintf(repmat('\b', 1, last_msg_len));
        fprintf('%s', msg);
        last_msg_len = length(msg);
    end
    
    if count >= total
        count = []; prev_percent = []; last_msg_len = [];
    end
end