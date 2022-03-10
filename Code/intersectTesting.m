if (0 <= py(k) && py(k) <= 0.4e-7)
    if (0.8e-7 <= px(k) && px(k) <= 1e-7)
        if (px_prev(k) < 0.8e-7)
            vx(k) = -vx(k);
            px(k) = 0.799e-7;
        end
        if (px_prev >= 0.8e-7)
            vy(k) = -vy(k);
            py(k) = 0.401e-7;
        end
    end
    if (1e-7 <= px(k) && px(k) <= 1.2e-7)
        if (px_prev(k) > 1.2e-7)
            vx(k) = -vx(k);
            px(k) = 1.201e-7;
        end
        if (px_prev(k) <= 1.2e-7)
            vy(k) = -vy(k);
            py(k) = 0.401e-7;
        end
    end
end

if (0.6e-7 <= py(k) && py(k) <= 1e-7)
    if (0.8e-7 <= px(k) && px(k) <= 1e-7)
        if (px_prev(k) < 0.8e-7)
            vx(k) = -vx(k);
            px(k) = 0.799e-7;
        end
        if (px_prev >= 0.8e-7)
            vy(k) = -vy(k);
            py(k) = 0.599e-7;
        end
    end
    if (1e-7 <= px(k) && px(k) <= 1.2e-7)
        if (px_prev(k) > 1.2e-7)
            vx(k) = -vx(k);
            px(k) = 1.201e-7;
        end
        if (px_prev(k) <= 1.2e-7)
            vy(k) = -vy(k);
            py(k) = 0.599e-7;
        end
    end
end
