function mtk_misc_print_to_tikz(x, y)
    for i=1:length(x)
        fprintf('(%f,%f)', x(i), y(i));
    end
    fprintf('\n');
end