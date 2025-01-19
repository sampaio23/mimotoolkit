function r = mtk_util_quantize(quantizer, y, factor)
    switch quantizer
        case '1bit'
            r = factor*(sign(real(y)) + 1i*sign(imag(y)));
        case 'unquant'
            r = y;
        otherwise
            error("Quantizer not implemented");
    end
end
