origin(rec::HyperRectangle) = rec.origin
widths(rec::HyperRectangle) = rec.widths
center(rec::HyperRectangle) = origin(rec) .+ widths(rec)
