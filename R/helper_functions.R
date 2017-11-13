convert_to_column_indices <- function(columns, block) {
    if (is.numeric(columns))
    {
        if (any(columns > NCOL(block)))
            warning("Some column indices exceed the number of columns ", 
                    "and were ignored.")
        return(columns[columns <= NCOL(block)])
    }
    # else
    indices <- match(columns, colnames(block))
    if (any(is.na(indices)))
        warning("Some column names could not be matched and were ignored.")
    return(indices[is.finite(indices)])
}