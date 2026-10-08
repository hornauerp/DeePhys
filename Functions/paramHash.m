function h = paramHash(params_struct)
% PARAMHASH  Deterministic MD5 hash of a parameter struct for cache keying.
%
%   h = paramHash(params_struct)
%
% Returns a 32-character hex string. Identical structs always produce the
% same hash, enabling parameter-keyed feature caching.
%
% The MD5 is computed in MATLAB (md5hex, ~0.1 ms) rather than through Java's
% MessageDigest: MATLAB R2025a+ ships without a JVM, and usejava('jvm') alone costs
% ~50 ms per call without one. The digest is identical to the former Java version, so
% existing RecordingIDs / UnitIDs are unchanged.

    h = md5hex(getByteStreamFromArray(params_struct));
end
