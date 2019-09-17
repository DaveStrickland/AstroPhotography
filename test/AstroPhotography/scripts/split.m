% @file split.m
% @brief Used to extract "truth" values for testing RawConv.split()
%
% @history 2019-09-14 dks : Initial version coded.
function split(r_file, g1_file, b_file, g2_file)
    format compact;
    more off;
    
    # Array size should be multiple of two to get full bayer mask in each image
    # otherwise we will distort any white-balance calculations.
    minrow = 451;
    maxrow = 464;
    mincol = 2851;
    maxcol = 2864;
    
    # index for sumval in out_stats.
    idx_sum = 8;
    
    [r_im,  r_stats]  = read_and_plot(r_file,  minrow, maxrow, mincol, maxcol);
    [g1_im, g1_stats] = read_and_plot(g1_file, minrow, maxrow, mincol, maxcol);
    [b_im,  b_stats]  = read_and_plot(b_file,  minrow, maxrow, mincol, maxcol);
    [g2_im, g2_stats] = read_and_plot(g2_file, minrow, maxrow, mincol, maxcol);
    
    # Calculate region-based white balance.
    # wb_in is sum from a given band, wb_out is multiplies for each band.
    wb_in = [r_stats(idx_sum) g1_stats(idx_sum) b_stats(idx_sum) g2_stats(idx_sum)]
    
    function [out_arr, out_stats] = read_and_plot(fname, minrow, maxrow, mincol, maxcol)
        % @brief Helper function to process a single image.
        %
        % @param minrow Mininum row number to include of image, inclusive.
        % @param maxrow Maximum row number to include of image, inclusive.
        % @param mincol Mininum column number to include of image, inclusive.
        % @param maxcol Maximum column number to include of image, inclusive.
        imdata = imread(fname);
        
        subplot(2,1,1);
        imagesc(imdata);
        title(fname);
        colorbar();
        axis image;
        
        out_arr = imdata(minrow:maxrow, mincol:maxcol);
        disp(sprintf('# %s %d %d %d %d', fname, ...
            minrow, maxrow, mincol, maxcol));
        out_stats = img_stats(out_arr);
        
        subplot(2,1,2);
        imagesc(out_arr);
        title(fname);
        colorbar();
        axis image;
        
        % Note an octave range a:b (1-based, inclusive) maps
        % to numpy (0-based, exclusive end) as a-1:b
        as_python(out_arr, '%d');
    end
    
    function as_python(arr, fmt)
        % Prints an octave array as a python list-of-lists,
        % suitable for use in python and/or numpy.
        %
        % @param arr Input array
        % @param fmt Printf-style numeric format, e.g. '%d' or '%.3f'
        if (ndims(arr) ~= 2)
            error(sprintf('as_python cannot handle %d-dimensional arrays.', ndims(arr)));
        end
        line_fmt = [fmt ', '];
        nrows = rows(arr);
        line_start_str=' [[';
        line_end_str='],';
        for ir = 1:nrows
            if (ir == 2)
                % Start of line
                line_start_str='  [';
            end
            if (ir == nrows)
                line_end_str=']]';
            end
        
            % Get a single row with commas
            line_str = num2str(arr(ir,:), line_fmt);
            
            % Trim off trailing whitespace and commas.
            line_str = strtrim(line_str);
            line_str = strtrunc(line_str, length(line_str)-1);
            
            % Add row braces. 
            line_str = [line_start_str line_str line_end_str];
            printf('%s\n', line_str);
        end
        % trailing end-of-array brace
    end
    
    function out_stats = img_stats(input_arr)
        % @brief Calculates and prints some simple statistics.
        nrows   = rows(input_arr);
        ncols   = columns(input_arr);
        nvals   = numel(input_arr);
        minval  = min(input_arr(:));
        maxval  = max(input_arr(:));
        sumval  = sum(input_arr(:));
        meanval = mean(input_arr(:));
        stdval  = std(input_arr(:));
        disp(sprintf('Array %dx%d has %d pixels', nrows, ncols, nvals));
        disp(sprintf('  meanval=%.f', meanval));
        disp(sprintf('  stdval=%.f', stdval));
        disp(sprintf('  minval=%.f', minval));
        disp(sprintf('  maxval=%.f', maxval));
        disp(sprintf('  sumval=%.f', sumval));
        out_stats = [nrows, ncols, nvals, minval, maxval, meanval, stdval, sumval];
    end
    
end
