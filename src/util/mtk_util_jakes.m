function eta = mtk_util_jakes(v, f, t)
    c = 3e8;
    eta = besselj(0, 2*pi*v*f*t/c);
end

