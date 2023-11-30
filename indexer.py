
#build the vector store
import chromadb
from llama_index import VectorStoreIndex, SimpleDirectoryReader
from llama_index.vector_stores import ChromaVectorStore
from llama_index import StorageContext

#ingest the data from downloaded articles
download_directory = '/Users/tscheidt/Documents/code/jdmgpt/downloaded_articles'
vector_store = ChromaVectorStore()
index = VectorStoreIndex(vector_store=vector_store)
reader = SimpleDirectoryReader(download_directory)
index.index(reader)

documents = SimpleDirectoryReader(download_directory).load_data()
index = VectorStoreIndex.from_documents(
    documents, service_context=service_context
)


#build persistent chroma index
chroma_client = chromadb.PersistentClient()
chroma_collection = chroma_client.create_collection("PCM_research")
vector_store = ChromaVectorStore(chroma_collection=chroma_collection)
storage_context = StorageContext.from_defaults(vector_store=vector_store)

# test query 
query_engine = index.as_query_engine()
response = query_engine.query("What did the author do growing up?")
print(response)