import random
import time
from typing import List, Optional
from Bio import Entrez


def ncbi_protein_id_extract(
        *,
        ncbi_email: str,
        ncbi_api_key: Optional[str],
        query: str,
) -> List[str]:
    """
        Extract NCBI protein IDs from an Entrez search query (batched).

        Args:
            ncbi_email: Required by NCBI usage policies / best practice.
            ncbi_api_key: Optional NCBI API key. If provided, higher request rate allowed.
            query: Entrez search term.

        Returns:
            List[str]: Protein IDs.
        """

    if not ncbi_email:
        raise ValueError(
            "No NCBI email credential found. Please add .env file with NCBI_EMAIL. NCBI_API_KEY is optional.")
    else:
        print("NCBI credentials found.")

    Entrez.email = ncbi_email
    Entrez.api_key = ncbi_api_key

    ret_max = 1000
    ret_start = 0

    # First esearch: save dataset in History Server (run snapshot)

    handle = Entrez.esearch(
        db='protein',
        term=query,
        retmax=0,
        idtype='acc',
        usehistory='y',
    )
    page = Entrez.read(handle)
    handle.close()

    protein_id_count = int(page['Count'])
    print(f'Found {protein_id_count} protein IDs.')

    webenv = page['WebEnv']
    query_key = page['QueryKey']

    print(f'History session: WebEnv={webenv}, QueryKey={query_key}.')

    protein_id_list: List[str] = []

    request_delay = 0.1 if ncbi_api_key else 0.34
    max_retries = 5

    # Consistent pagination with WebEnv and query_key
    while ret_start < protein_id_count:

        for attempt in range(max_retries):

            try:
                handle = Entrez.esearch(
                    db='protein',
                    term=query,
                    idtype='acc',
                    usehistory='y',
                    retmax=ret_max,
                    retstart=ret_start,
                    WebEnv=webenv,
                    query_key=query_key,
                )
                page = Entrez.read(handle)
                handle.close()

                protein_ids_this_batch = page['IdList']
                protein_id_list.extend(protein_ids_this_batch)

                ret_start += ret_max

                if protein_id_count > 0:
                    progress = min(ret_start / protein_id_count, 1.0)
                    print(f'Extracted {len(protein_id_list)} protein IDs ({progress:.2%}).')

                jitter = random.uniform(0.01, 0.05)
                time.sleep(request_delay + jitter)

                break

            except Exception as e:

                if attempt == max_retries - 1:
                    raise RuntimeError(
                        f"Failed to extract protein IDs after {max_retries} attempts at retstart={ret_start}: {e}"
                    )

                backoff = (2 ** attempt) * request_delay
                print(f"Retry {attempt + 1}/{max_retries} after error: {e}")
                time.sleep(backoff)

    print(f"Extracted {len(protein_id_list)} protein IDs.")
    return protein_id_list
