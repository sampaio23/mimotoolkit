function x_hat = mtk_util_demodulate(modulation, x_tilde)
    switch modulation
        case 'qpsk'
            x_hat = 1/sqrt(2)*(sign(real(x_tilde)) + 1i*sign(imag(x_tilde)));
        otherwise
            error('Modulation not implemented');
    end
end

