function t_str = formatTime(seconds)
    % Formats seconds into HH:MM:SS or MM:SS
    if seconds < 0, seconds = 0; end
    h = floor(seconds / 3600);
    m = floor(mod(seconds, 3600) / 60);
    s = floor(mod(seconds, 60));
    
    if h > 0
        t_str = sprintf('%02dh %02dm %02ds', h, m, s);
    else
        t_str = sprintf('%02dm %02ds', m, s);
    end
end