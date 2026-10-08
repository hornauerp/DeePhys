function h = md5hex(bytes)
% MD5HEX  RFC 1321 MD5 of a uint8 vector, as a lower-case hex string.
%   Pure MATLAB (no JVM). Bit operations run on uint32; sums are taken in double and
%   reduced mod 2^32 (exact, all terms < 2^53).
    msg = uint8(bytes(:));
    n   = numel(msg);
    len_bytes = typecast(uint64(n) * 8, 'uint8');          % little-endian bit length
    msg = [msg; uint8(128); zeros(mod(55 - n, 64), 1, 'uint8'); len_bytes(:)];
    words = reshape(typecast(msg, 'uint32'), 16, []);

    s = [7 12 17 22 7 12 17 22 7 12 17 22 7 12 17 22, ...
         5  9 14 20 5  9 14 20 5  9 14 20 5  9 14 20, ...
         4 11 16 23 4 11 16 23 4 11 16 23 4 11 16 23, ...
         6 10 15 21 6 10 15 21 6 10 15 21 6 10 15 21];
    K = floor(abs(sin(1:64)) * 2^32);
    g_idx = [0:15, mod(5 * (16:31) + 1, 16), mod(3 * (32:47) + 5, 16), mod(7 * (48:63), 16)] + 1;
    M32 = 2^32;

    a0 = uint32(1732584193); b0 = uint32(4023233417); c0 = uint32(2562383102); d0 = uint32(271733878);
    for chunk = 1:size(words, 2)
        Mw = double(words(:, chunk));
        A = a0; B = b0; C = c0; D = d0;
        for i = 1:64
            if i <= 16
                F = bitor(bitand(B, C), bitand(bitcmp(B), D));
            elseif i <= 32
                F = bitor(bitand(D, B), bitand(bitcmp(D), C));
            elseif i <= 48
                F = bitxor(bitxor(B, C), D);
            else
                F = bitxor(C, bitor(B, bitcmp(D)));
            end
            F = uint32(mod(double(F) + double(A) + K(i) + Mw(g_idx(i)), M32));
            A = D;  D = C;  C = B;
            rot = bitor(bitshift(F, s(i)), bitshift(F, s(i) - 32));   % rotate left
            B = uint32(mod(double(B) + double(rot), M32));
        end
        a0 = uint32(mod(double(a0) + double(A), M32)); b0 = uint32(mod(double(b0) + double(B), M32));
        c0 = uint32(mod(double(c0) + double(C), M32)); d0 = uint32(mod(double(d0) + double(D), M32));
    end
    h = sprintf('%02x', typecast([a0 b0 c0 d0], 'uint8'));
end
