function H_real = mtk_util_mat_real(H)
    H_real = [real(H) -imag(H); imag(H) real(H)];
end


    
