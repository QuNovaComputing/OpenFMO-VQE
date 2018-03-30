
#include <string.h>
#include <falanx.h>

const int const g_ofmo_storage_data_ids_null[6] = { -1, -1, -1, -1, -1, -1 };

/*
 * static変数
 */
static falanx_context_t* CONTEXT;
static int DATA_IDS[6];


/*
 * falanx_context_t* を保存する.
 */
void
ofmo_storage_set_context(falanx_context_t* context)
{
    CONTEXT = context;
}

falanx_context_t*
ofmo_storage_get_context(void)
{
    return CONTEXT;
}

void
ofmo_storage_set_data_ids(const int data_ids[6])
{
    memcpy(DATA_IDS, data_ids, sizeof(int)*6);
}

void
ofmo_storage_get_data_ids(int data_ids[6])
{
    memcpy(data_ids, DATA_IDS, sizeof(int)*6);
}


