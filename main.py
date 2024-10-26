import pandas as pd


def read_file(path, sep=",", **kwargs):
    """
    Função para ler arquivos com diagnóstico.
    """
    try:
        df = pd.read_csv(path, sep=sep, dtype=str, **kwargs)
        print(f"Arquivo '{path}' lido com sucesso. Número de linhas: {len(df)}")
        print(f"Colunas em '{path}': {df.columns.tolist()}\n")
        duplicados = df.duplicated().sum()
        if duplicados > 0:
            df = df.drop_duplicates(keep=False, inplace=True)
            print(f"Linhas duplicadas removidas: {duplicados}")
        else:
            print("Nenhuma linha duplicada encontrada.")

        return df

    except pd.errors.ParserError as e:
        print(f"Erro ao ler o arquivo '{path}': {e}")
        df = pd.read_csv(path, sep=sep, dtype=str, on_bad_lines="skip", **kwargs)
        print(
            f"Arquivo '{path}' lido com sucesso após ignorar linhas com erros. Número de linhas: {len(df)}"
        )
        print(f"Colunas em '{path}': {df.columns.tolist()}\n")

        duplicados = df.duplicated().sum()
        if duplicados > 0:
            df = df.drop_duplicates(keep=False, inplace=True)
            print(f"Linhas duplicadas removidas: {duplicados}")
        else:
            print("Nenhuma linha duplicada encontrada.")

        return df


def transform_peptide(peptide):
    """
    Função para converter a sequência de peptídeo para maiúsculas
    e substituir 'I' e 'L' por 'X'.
    """
    if pd.isnull(peptide):
        return peptide
    return peptide.upper().replace("I", "X").replace("L", "X")


conversion_table = read_file("./processed_files/tabela_np_gi.csv")
conversion_table.set_index(["NP"], inplace=True)


def np_to_gi(np):
    return conversion_table.loc[np]["GI"]


evidence_pepvar_file = "./processed_files/evidence_dbpepvar.csv"
evidence_refseq_file = "./processed_files/evidence_refseq.csv"
missense_file = "./processed_files/missense.csv"


def main():
    evidence_pepvar = read_file(evidence_pepvar_file)
    evidence_refseq = read_file(evidence_refseq_file)
    missense = read_file(missense_file)

    evidence_refseq["Sequence"] = evidence_refseq["Sequence"].apply(transform_peptide)
    evidence_pepvar["Sequence"] = evidence_pepvar["Sequence"].apply(transform_peptide)

    missense["PepRef"] = missense["PepRef"].apply(transform_peptide)
    missense["PepMut"] = missense["PepMut"].apply(transform_peptide)

    print(evidence_refseq.columns.tolist())
    print(evidence_pepvar.columns.tolist())
    evidence_pepvar["Leading Razor Protein"] = evidence_pepvar[
        "Leading Razor Protein"
    ].apply(np_to_gi)

    pd.merge(evidence_refseq, evidence_pepvar, left_on=["Raw file", "MS/MS Scan Number"], right_on=['Raw file', 'MS MS Scan Number'])

    print("Fazendo merge dos evidences")
    merged_df = pd.merge(
        evidence_pepvar,
        evidence_refseq,
        on=["Raw file", "Leading Razor Protein"],
        how="inner",
        suffixes=("_dbPepVar", "_refSeq"),
    )

    print("Fazendo merge do missense")
    output = pd.merge(
        merged_df,
        missense,
        left_on=["Sequence_dbPepVar", "Sequence_refSeq"],
        right_on=["PepMut", "PepRef"],
        how="inner",
    )

    print(output)

    output.to_csv("output.csv", index=False)


if __name__ == "__main__":
    main()
