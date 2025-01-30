function H = mtk_util_mat_complex(H_real)
    rows = size(H_real, 1);
    cols = size(H_real, 2);
    H = H_real(1 : rows/2, 1 : cols/2) + 1i * H_real(rows/2 + 1 : rows, 1 : cols/2);
end


    
