import streamlit as st
from llama_index import VectorStoreIndex, ServiceContext, Document
from llama_index.llms import OpenAI
import openai
import os
import chromadb
from llama_index.vector_stores import ChromaVectorStore
from llama_index import SimpleDirectoryReader, StorageContext

st.set_page_config(page_title="Chat with open access JDM research", page_icon="🦙", layout="centered", initial_sidebar_state="auto", menu_items=None)
openai.api_key = os.environ["OPENAI_API_KEY"]
st.title("Chat with the JDM research 💬🦙")
st.info("Say something here", icon="📃")
         
if "messages" not in st.session_state.keys(): # Initialize the chat messages history
    st.session_state.messages = [
        {"role": "assistant", "content": "Ask me a question about JDM research!"}
    ]

@st.cache_resource(show_spinner=False)

def load_data():
    with st.spinner(text="Loading and indexing the JDM research docs – hang tight! This should take 1-2 minutes."):
        # build persistent chroma client
        chroma_client = chromadb.PersistentClient(path="/Users/tscheidt/Documents/code/jdmgpt/chroma_db")
        chroma_collection = chroma_client.get_or_create_collection("gptbot") # create a new collection if doesn't exist
        print(chroma_collection.name)
        print(chroma_collection.count())

        reader = SimpleDirectoryReader(input_dir="./sample_articles", recursive=True)
        docs = reader.load_data()
        #print(docs[0])

        vector_store = ChromaVectorStore(chroma_collection=chroma_collection)
        storage_context = StorageContext.from_defaults(vector_store=vector_store)
        service_context = ServiceContext.from_defaults(llm=OpenAI(model="gpt-3.5-turbo", temperature=0.5, system_prompt="You are an expert medical researcher autoimmune conditions with a specialization in juvenile dermatomyositis (JDM) and your job is to answer questions from medical practitioners and patients to help understand the causes and potential treatments of JDM. Assume that JDM means juvenile dermatomyositis and all questions are related to the JDM Research Library. Keep your answers technical and based on facts – do not hallucinate features."))
        
        index = VectorStoreIndex.from_documents(docs, storage_context=storage_context, service_context=service_context)
        #print(index[0:10])
        return index

index = load_data()

# chat engine initialization - details of options here: https://docs.llamaindex.ai/en/latest/module_guides/deploying/chat_engines/usage_pattern.html#
if "chat_engine" not in st.session_state.keys(): # Initialize the chat engine
        st.session_state.chat_engine = index.as_chat_engine(chat_mode="openai", verbose=True)

if prompt := st.chat_input("Your question"): # Prompt for user input and save to chat history
    st.session_state.messages.append({"role": "user", "content": prompt})

for message in st.session_state.messages: # Display the prior chat messages
    with st.chat_message(message["role"]):
        st.write(message["content"])

# If last message is not from assistant, generate a new response
if st.session_state.messages[-1]["role"] != "assistant":
    with st.chat_message("assistant"):
        with st.spinner("Thinking..."):
            response = st.session_state.chat_engine.chat(prompt)
            st.write(response.response)
            message = {"role": "assistant", "content": response.response}
            st.session_state.messages.append(message) # Add response to message history