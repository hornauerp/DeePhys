function h = paramHash(params_struct)
% PARAMHASH  Deterministic MD5 hash of a parameter struct for cache keying.
%
%   h = paramHash(params_struct)
%
% Returns a 32-character hex string. Identical structs always produce the
% same hash, enabling parameter-keyed feature caching.

    bytes = getByteStreamFromArray(params_struct);
    md = java.security.MessageDigest.getInstance('MD5');
    md.update(bytes);
    hash_bytes = typecast(md.digest(), 'uint8');
    h = sprintf('%02x', hash_bytes);
end
