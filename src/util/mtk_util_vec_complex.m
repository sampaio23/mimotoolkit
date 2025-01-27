function x = mtk_util_vec_complex(x_real)
    rows = size(x_real, 1);
    x = x_real(1 : rows/2, :) + 1i * x_real(rows/2 + 1 : rows, :);
end


    
