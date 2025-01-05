function r = mtk_util_quantize(y, factor)
    r = factor*(sign(real(y)) + 1i*sign(imag(y)));
end
