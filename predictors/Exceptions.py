from constants import included_genes


class InvalidFormatError(Exception):
     def __init__(self, message):
         self.message = message
         super(InvalidFormatError, self).__init__(self.message)


class VariantFormatException(Exception):
    def __init__(self, variant):
        self.variant = variant
        self.message = f"Invalid Variant Format for the variant : {variant}"
        super(VariantFormatException, self).__init__(self.message)

class VariantNotInGeneException(Exception):
    def __init__(self, variant, chrom, start, end, gene):
        self.message = f"""
        Variant not in bounds of specified gene
        Variant {variant} is not in {gene}
        {gene} Assembly 38 coordiantes are: {chrom}:g.{start} - {chrom}:g.{end}
        """
        super(VariantNotInGeneException, self).__init__(self.message)

class InvalidGeneException(Exception):
    def __init__(self, input_gene):
        gene_list = ", ".join(included_genes)
        self.message = f"""
         {input_gene} predictive model not available. Available gene models are :
         {gene_list}
        """
        super(InvalidGeneException, self).__init__(self.message)

class NotInitializedException(Exception):
    def __init__(self):
        self.message = f"""
         Model not initialized
        """
        super(NotInitializedException, self).__init__(self.message)

class ReferenceSequenceException(Exception):
    def __init__(self, correct_allele, observed_allele, variant):
        self.message = f"""
            Refernce allele {correct_allele} does not match the reported allele {observed_allele}
            in the variant {variant}.
        """
        super(ReferenceSequenceException, self).__init__(self.message)
