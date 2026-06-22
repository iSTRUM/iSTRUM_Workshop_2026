function idx = find_index_where(x, val)
    % find the index of the value in x that is closest to val
    [~, idx] = min(abs(x - val));
end